import json
import sys
import csv
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
from rocketcea.cea_obj import CEA_Obj, add_new_fuel

# n-Butanol (1-butanol) card
# h,cal = heat of formation = -159,600 cal/mol (liquid, 298.15 K)
# rho = 0.810 g/cc
card_str = """
fuel  C4H10O(L)  C 4  H 10  O 1
      wt%=100.00  h,cal=-78196.0  t(k)=298.15  rho=0.810
"""

add_new_fuel('nButanol', card_str)
# -----------------------------
# Physics Helpers
# -----------------------------
INtoM = 0.0254
def K_to_F(Tk: float) -> float: return (Tk - 273.15) * 9.0 / 5.0 + 32.0
def R_to_K(TR: float) -> float: return TR * (5.0 / 9.0)
def psi_to_Pa(psi: float) -> float: return psi * 6894.757293168

def thomas_solve(a, d, b, c):
    n_ = len(d)
    dp, cpv = d.astype(float).copy(), c.astype(float).copy()
    for i in range(1, n_):
        m = b[i - 1] / dp[i - 1]
        dp[i] -= m * a[i - 1]
        cpv[i] -= m * cpv[i - 1]
    x = np.empty(n_, dtype=float)
    x[-1] = cpv[-1] / dp[-1]
    for i in range(n_ - 2, -1, -1):
        x[i] = (cpv[i] - a[i] * x[i + 1]) / dp[i]
    return x

# -----------------------------
# CEA Station Property Lookup
# -----------------------------
VALID_STATIONS = {"chamber_section", "throat_section"}

def get_cea_station_properties(cea: CEA_Obj, cea_cfg: dict, station: str) -> dict:
    """
    Retrieve gas temperature, transport properties, and geometry reference
    from RocketCEA for the specified station.

    Supported stations
    ------------------
    "chamber_section"
        - T_g  : chamber (combustion) temperature  → get_Temperatures()[0]
        - gamma: chamber gamma                      → get_Chamber_MolWt_gamma()
        - transport (Cp, visc, cond, Pr)            → get_Chamber_Transport()
        - Ma   : defaults to 0.3 (subsonic chamber)

    "throat_section"
        - T_g  : throat static temperature         → get_Temperatures()[1]
        - gamma: throat gamma                       → get_Throat_MolWt_gamma()
        - transport (Cp, visc, cond, Pr)            → get_Throat_Transport()
        - Ma   : 1.0 (sonic throat)

    Returns a dict with keys:
        T_g, gam, cp_g, mu, Pr, Ma_ref
    """
    if station not in VALID_STATIONS:
        raise ValueError(
            f"Unknown station '{station}'. "
            f"Must be one of: {sorted(VALID_STATIONS)}"
        )

    Pc = cea_cfg["Pc_chamber_psia"]
    MR = cea_cfg["MR"]

    # --- Temperature (index 0 = chamber, 1 = throat, 2 = exit) ---
    all_temps = cea.get_Temperatures(Pc=Pc, MR=MR, frozen=False)

    if station == "chamber_section":
        T_g_R  = float(all_temps[0])
        _, gam = cea.get_Chamber_MolWt_gamma(Pc=Pc, MR=MR)
        Cp_c, visc_c, cond_c, Pr = cea.get_Chamber_Transport(Pc=Pc, MR=MR)
        Ma_ref = cea_cfg.get("Ma_chamber", 0.3)

    else:  # throat_section
        T_g_R  = float(all_temps[1])
        _, gam = cea.get_Throat_MolWt_gamma(Pc=Pc, MR=MR)
        Cp_c, visc_c, cond_c, Pr = cea.get_Throat_Transport(Pc=Pc, MR=MR)
        Ma_ref = cea_cfg.get("Ma_throat", 1.0)

    T_g  = R_to_K(T_g_R)
    cp_g = float(Cp_c) * 4184.0   # cal/(g·K) → J/(kg·K)
    mu   = float(visc_c) * 1e-4   # milliPoise → Pa·s

    return dict(T_g=T_g, gam=float(gam), cp_g=cp_g, mu=mu, Pr=float(Pr), Ma_ref=Ma_ref)

# -----------------------------
# Simulation Management
# -----------------------------
class SimulationResult:
    def __init__(self, name, x_in, t_melt_f, frames_data, times, melt_time, edge_log, pc_psia, station):
        self.name, self.x_in, self.t_melt_f = name, x_in, t_melt_f
        self.frames_data, self.times, self.melt_time, self.edge_log = frames_data, times, melt_time, edge_log
        self.pc_psia = pc_psia
        self.station = station

def load_json_safe(path: Path):
    try:
        with open(path, "r") as f:
            return json.load(f)
    except Exception as e:
        raise SystemExit(f"Error loading {path.name}: {e}")

def run_simulation(run_entry: dict, setup_path: Path):
    config     = load_json_safe(setup_path / "config.json")
    materials  = load_json_safe(setup_path / "materials.json")
    geometries = load_json_safe(setup_path / "geometry.json")

    mat  = materials[run_entry["material_id"]]
    geom = geometries[run_entry["geometry_id"]]
    cea_cfg, solv_cfg, env_cfg = config["cea"], config["solver"], config["environment"]

    # Station is fully derived from geometry_id
    station = run_entry["geometry_id"]

    L, n = float(geom["L_in"]) * INtoM, int(solv_cfg["n"])
    k_w, rho, cp_w, T_melt = float(mat["k_w"]), float(mat["rho"]), float(mat["cp_wall"]), float(mat["T_melt_K"])
    h_nat, T_amb, T_init = float(env_cfg["h_nat"]), float(env_cfg["T_amb_K"]), float(env_cfg["T_init_K"])
    tf, dt = float(solv_cfg["tf"]), float(solv_cfg["dt"])

    Thr_D_in     = float(geom["Throat_D_in"])
    section_D_in = float(geom["D_in"])
    A            = np.pi * (section_D_in * INtoM / 2) ** 2
    Thr_A        = np.pi * (Thr_D_in    * INtoM / 2) ** 2

    # --- CEA setup & station-specific properties ---
    cea = CEA_Obj(oxName=cea_cfg["oxName"], fuelName=cea_cfg["fuelName"])

    props = get_cea_station_properties(cea, cea_cfg, station)
    T_g  = props["T_g"]
    gam  = props["gam"]
    cp_g = props["cp_g"]
    mu   = props["mu"]
    Pr   = props["Pr"]
    Ma   = props["Ma_ref"]

    print(
        f"  Station : {station}\n"
        f"  T_g     : {T_g:.1f} K  ({K_to_F(T_g):.1f} °F)\n"
        f"  gamma   : {gam:.4f}\n"
        f"  cp_g    : {cp_g:.1f} J/(kg·K)\n"
        f"  mu      : {mu:.3e} Pa·s\n"
        f"  Pr      : {Pr:.4f}\n"
        f"  Ma_ref  : {Ma:.2f}"
    )

    Pc_Pa = psi_to_Pa(cea_cfg["Pc_chamber_psia"])
    Cstar = float(cea.get_Cstar(Pc=cea_cfg["Pc_chamber_psia"], MR=cea_cfg["MR"])) * 0.3048

    dy  = L / n
    phi = T_init * np.ones(n)

    # Bartz base heat-transfer coefficient (referenced to throat diameter)
    area_ratio_term = (Thr_A / A) ** 0.9 if station == "chamber_section" else 1.0

    hG_base = (
        (0.026 / ((Thr_D_in * INtoM) ** 0.2))
        * ((mu ** 0.2 * cp_g) / (Pr ** 0.6))
        * (Pc_Pa / Cstar) ** 0.8
        * area_ratio_term
    )

    # Adiabatic wall temperature at the chosen station
    Taw = T_g * (1.0 + Pr**(1/3) * (gam - 1.0)/2.0 * Ma**2)
    DD, ap0 = k_w / dy, rho * cp_w * dy / dt

    history, times, edge_log, melt_time = [], [], [], None

    for step in range(int(tf / dt) + 1):
        t     = step * dt
        Trat  = 1.0 + (gam - 1.0) / 2.0 * Ma ** 2
        sigma = 1.0 / (((0.5 * (phi[0] / T_g) * Trat + 0.5) ** 0.68) * (Trat ** 0.12))
        hG    = hG_base * sigma

        d      = (ap0 + 2.0 * DD) * np.ones(n)
        d[0]   = ap0 / 2 + DD + hG
        d[-1]  = ap0 / 2 + DD + h_nat
        a      = -DD * np.ones(n - 1)
        b      = -DD * np.ones(n - 1)
        c      = ap0 * phi
        c[0]   = (ap0 / 2) * phi[0] + hG * Taw
        c[-1]  = (ap0 / 2) * phi[-1] + h_nat * T_amb
        phi    = thomas_solve(a, d, b, c)

        if melt_time is None and phi[0] >= T_melt:
            melt_time = t
        if step % max(1, int(round((1 / 10) / dt))) == 0:
            history.append(phi.copy())
            times.append(t)
            edge_log.append([t, K_to_F(phi[0])])

    return SimulationResult(
        run_entry["run_name"],
        np.linspace(0, L, n) / INtoM,
        K_to_F(T_melt),
        history, times, melt_time, edge_log,
        cea_cfg["Pc_chamber_psia"],
        station,
    )

def render_gif(data_list, filename, is_combined=False):
    fig, ax = plt.subplots(figsize=(10, 6))

    max_wall_thickness = max(r.x_in[-1] for r in data_list)
    ax.set_xlim(0.0, max_wall_thickness)
    ax.set_ylim(K_to_F(280), 3500)
    ax.set_xlabel("Depth (in)")
    ax.set_ylabel("Temperature (°F)")
    ax.grid(True)

    lines = []
    for r in data_list:
        lbl = f"{r.name} [{r.station}]"
        if r.melt_time:
            lbl += f" (Melt: {r.melt_time:.2f}s)"
        l, = ax.plot(r.x_in, K_to_F(r.frames_data[0]), label=lbl, linewidth=2)
        ax.axhline(y=r.t_melt_f, color=l.get_color(), linestyle='--', alpha=0.4)
        lines.append(l)

    ax.legend(loc='upper right', fontsize='small')
    frames = []

    max_frames = max(len(r.frames_data) for r in data_list)
    for i in range(max_frames):
        for idx, r in enumerate(data_list):
            f_idx = min(i, len(r.frames_data) - 1)
            lines[idx].set_ydata(K_to_F(r.frames_data[f_idx]))

        curr_t = data_list[0].times[min(i, len(data_list[0].times) - 1)]

        if is_combined:
            ax.set_title(f"Combined Analysis | Time: {curr_t:.2f} s")
        else:
            res = data_list[0]
            ax.set_title(
                f"{res.name} [{res.station}] | Time: {curr_t:.2f} s | "
                f"Pc: {res.pc_psia:.0f} psia | "
                f"Wall: {K_to_F(res.frames_data[min(i, len(res.frames_data)-1)][0]):.1f} °F"
            )

        fig.canvas.draw()
        frames.append(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy())

    imageio.mimsave(filename, frames, fps=30)
    plt.close(fig)

def main():
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python script.py <Project_Folder>")
    project_dir = Path(sys.argv[1]).expanduser().resolve()
    setup_dir   = project_dir / "Setup"
    run_list    = load_json_safe(project_dir / "run.json")

    results_base = project_dir / f"{project_dir.name}_Results"
    results_base.mkdir(parents=True, exist_ok=True)

    mode = "individual"
    if len(run_list) > 1:
        if input(f"Run {len(run_list)} cases together? (t/i): ").lower() == 't':
            mode = "together"

    all_res = []
    for entry in run_list:
        print(f"\nSolving: {entry['run_name']}")
        res     = run_simulation(entry, setup_dir)
        out_dir = results_base / entry["run_name"]
        out_dir.mkdir(parents=True, exist_ok=True)

        with open(out_dir / f"{entry['run_name']}_log.csv", "w", newline="") as f:
            csv.writer(f).writerows([["t_s", "wall_surface_F"]] + res.edge_log)

        render_gif([res], out_dir / f"{entry['run_name']}.gif", is_combined=False)
        if mode == "together":
            all_res.append(res)

        if entry.get("expansion_flag", 0) == 1:
            mat   = load_json_safe(setup_dir / "materials.json")[entry["material_id"]]
            E     = float(mat["E_tension_Pa"])
            nu    = float(mat["poisson"])
            alpha = float(mat["CTE_per_K"])
            T_ref = float(mat["CTE_ref_temp_K"])

            T_final        = np.mean(res.frames_data[-1])
            thermal_strain = alpha * (T_final - T_ref)
            thermal_stress = E * thermal_strain
            print(f"  Thermal Strain : {thermal_strain:.6f}")
            print(f"  Thermal Stress : {thermal_stress / 1e6:.3f} MPa")

    if mode == "together" and all_res:
        print("\nCreating combined GIF...")
        render_gif(all_res, results_base / "Combined_Analysis.gif", is_combined=True)


if __name__ == "__main__":
    main()