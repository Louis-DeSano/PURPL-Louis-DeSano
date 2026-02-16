"""
1-D Convection Heat Transfer Solver (Python) — file-driven inputs

Reads:
- CEA inputs + geometry + solver settings from a JSON config file
- Wall / boundary properties from the same JSON (keeps everything in one place)

Usage:
  python heat_transfer_1d.py config.json

Outputs (relative to this script unless overridden):
  ./Results/Heat_Transfer_Simulation.gif
  ./Results/Edge_Temperature_Log.csv
"""

import json
import sys
import csv
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio
from rocketcea.cea_obj import CEA_Obj


# -----------------------------
# Helpers
# -----------------------------
INtoM = 0.0254

def K_to_F(Tk: float) -> float:
    return (Tk - 273.15) * 9.0 / 5.0 + 32.0

def R_to_K(TR: float) -> float:
    return TR * (5.0 / 9.0)

def psi_to_Pa(psi: float) -> float:
    return psi * 6894.757293168

def get_script_dir() -> Path:
    # Works in scripts and notebooks
    try:
        return Path(__file__).resolve().parent
    except NameError:
        return Path.cwd()

def load_config(config_path: Path) -> dict:
    with open(config_path, "r") as f:
        cfg = json.load(f)

    # Minimal validation
    if "cea" not in cfg:
        raise ValueError("Config missing required top-level key: 'cea'")
    if "wall" not in cfg:
        raise ValueError("Config missing required top-level key: 'wall'")
    if "geometry" not in cfg:
        raise ValueError("Config missing required top-level key: 'geometry'")
    if "solver" not in cfg:
        raise ValueError("Config missing required top-level key: 'solver'")
    return cfg

def iter_config_paths(arg: str):
    p = Path(arg).expanduser().resolve()

    if p.is_file():
        if p.suffix.lower() != ".json":
            raise SystemExit(f"Expected a .json file, got: {p.name}")
        return [p]

    if p.is_dir():
        cfgs = sorted(p.glob("*.json"))
        if not cfgs:
            raise SystemExit(f"No .json files found in folder: {p}")
        return cfgs

    raise SystemExit(f"Path not found: {p}")


# -----------------------------
# Thomas solver
# -----------------------------
def thomas_solve(a, d, b, c):
    n_ = len(d)
    dp = d.astype(float).copy()
    cpv = c.astype(float).copy()

    for i in range(1, n_):
        m = b[i - 1] / dp[i - 1]
        dp[i] -= m * a[i - 1]
        cpv[i] -= m * cpv[i - 1]

    x = np.empty(n_, dtype=float)
    x[-1] = cpv[-1] / dp[-1]
    for i in range(n_ - 2, -1, -1):
        x[i] = (cpv[i] - a[i] * x[i + 1]) / dp[i]
    return x


def run_one_config(config_path: Path):
    cfg = load_config(config_path)

    # -----------------------------
    # Unpack config
    # -----------------------------
    cea_cfg = cfg["cea"]
    wall_cfg = cfg["wall"]
    geom_cfg = cfg["geometry"]
    solv_cfg = cfg["solver"]
    out_cfg = cfg.get("output", {})

    station_flag = str(cea_cfg.get("property_station", "throat")).lower()

    if station_flag not in ["throat", "chamber"]:
        raise ValueError("cea.property_station must be 'throat' or 'chamber'")

    USE_THROAT_PROPERTIES = (station_flag == "throat")

    Pc_chamber_psia = float(cea_cfg["Pc_chamber_psia"])
    MR = float(cea_cfg["MR"])
    oxName = str(cea_cfg["oxName"])
    fuelName = str(cea_cfg["fuelName"])
    eps = float(cea_cfg.get("eps", 1.0))
    frozen = int(cea_cfg.get("frozen", 0))  # keep your current default (0)

    # Geometry
    D = float(geom_cfg["D_in"]) * INtoM
    A = np.pi * (D / 2.0) ** 2
    Ma = float(geom_cfg.get("Ma", 1.0))

    Throat_D = float(geom_cfg.get("Throat_D_in", geom_cfg["D_in"])) * INtoM
    Throat_A = np.pi * (Throat_D / 2.0) ** 2

    # Wall / boundaries
    L = float(wall_cfg["L_in"]) * INtoM
    n = int(wall_cfg["n"])
    h_nat = float(wall_cfg["h_nat"])
    T_amb = float(wall_cfg["T_amb_K"])
    k_w = float(wall_cfg["k_w"])
    rho = float(wall_cfg["rho"])
    cp_wall = float(wall_cfg["cp_wall"])
    T_init = float(wall_cfg.get("T_init_K", T_amb))
    T_melt = float(wall_cfg.get("T_melt_K", 1673.0))

    # Solver
    tf = float(solv_cfg["tf"])
    dt = float(solv_cfg["dt"])

    # Output / logging / gif
    GIF_FPS = int(out_cfg.get("gif_fps", 30))
    PROGRESS_DT = float(out_cfg.get("progress_dt", 0.25))
    results_subdir = str(out_cfg.get("results_dir", "Results"))
    gif_name = str(out_cfg.get("gif_name", "Heat_Transfer_Simulation.gif"))
    csv_name = str(out_cfg.get("csv_name", "Edge_Temperature_Log.csv"))

# -----------------------------
# Results path (named after JSON)
# -----------------------------
    BASE_DIR = get_script_dir()

    config_stem = config_path.stem  # filename without .json

    RESULTS_DIR = (BASE_DIR / "Results" / config_stem).resolve()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    OUT_GIF = RESULTS_DIR / gif_name
    OUT_CSV = RESULTS_DIR / csv_name

    # -----------------------------
    # RocketCEA: pull properties (no try/except)
    # -----------------------------
    cea = CEA_Obj(oxName=oxName, fuelName=fuelName)

    temps_R = cea.get_Temperatures(Pc=Pc_chamber_psia, MR=MR, eps=eps, frozen=frozen)
    T_chamber_K = R_to_K(float(temps_R[0]))
    T_throat_K = R_to_K(float(temps_R[1]))

    _, gam_ch = cea.get_Chamber_MolWt_gamma(Pc=Pc_chamber_psia, MR=MR, eps=eps)
    _, gam_t = cea.get_Throat_MolWt_gamma(Pc=Pc_chamber_psia, MR=MR, eps=eps)
    gam_ch = float(gam_ch)
    gam_t = float(gam_t)

    Cp_ch, visc_ch, cond_ch, Pr_ch = cea.get_Chamber_Transport(
        Pc=Pc_chamber_psia, MR=MR, eps=eps, frozen=frozen
    )
    Cp_t, visc_t, cond_t, Pr_t = cea.get_Throat_Transport(
        Pc=Pc_chamber_psia, MR=MR, eps=eps, frozen=frozen
    )

    # Convert to SI
    cp_ch = float(Cp_ch) * 4184.0
    cp_t = float(Cp_t) * 4184.0
    mu_ch = float(visc_ch) * 1.0e-4
    mu_t = float(visc_t) * 1.0e-4
    k_ch = float(cond_ch) * 0.4184
    k_t = float(cond_t) * 0.4184
    Pr_ch = float(Pr_ch)
    Pr_t = float(Pr_t)

    # Bartz-required chamber Pc and chamber c*
    Pc_chamber_Pa = psi_to_Pa(Pc_chamber_psia)
    Cstar_chamber_m_s = float(cea.get_Cstar(Pc=Pc_chamber_psia, MR=MR)) * 0.3048

    # Toggle only transport inputs (Pc and c* remain chamber)
    if USE_THROAT_PROPERTIES:
        station = "THROAT_PROPERTIES"
        T_g = T_throat_K
        Pr = Pr_t
        mu = mu_t
        k_g = k_t
        cp_g = cp_t
        gam = gam_t
    else:
        station = "CHAMBER_PROPERTIES"
        T_g = T_chamber_K
        Pr = Pr_ch
        mu = mu_ch
        k_g = k_ch
        cp_g = cp_ch
        gam = gam_ch

    print(f"RocketCEA mode: {station}")
    print(f"  Bartz uses: Pc(chamber) = {Pc_chamber_psia:.3f} psia, c*(chamber) = {Cstar_chamber_m_s:.2f} m/s")
    print(f"  Transport uses: T_g={T_g:.2f} K, gamma={gam:.5f}, Pr={Pr:.5f}, mu={mu:.6e} Pa*s, cp={cp_g:.2f} J/(kg*K)")
    print()

    # -----------------------------
    # Discretization
    # -----------------------------
    dy = L / n
    phi = T_init * np.ones(n, dtype=float)

    # -----------------------------
    # Bartz convection coefficient (STANDARD)
    # -----------------------------
    hG_base = (
        (0.026 / (Throat_D ** 0.2))
        * ((mu ** 0.2 * cp_g) / (Pr ** 0.6))
        * (Pc_chamber_Pa / Cstar_chamber_m_s) ** 0.8
        * (Throat_A / A) ** 0.9
    )

    Taw = T_g * (1.0 + (gam - 1.0) / 2.0 * (Ma ** 2) * (Pr ** (1.0 / 3.0)))

    DD = k_w / dy
    ap0 = rho * cp_wall * dy / dt

    # -----------------------------
    # GIF / progress controls
    # -----------------------------
    FRAME_STRIDE = int(round((1.0 / GIF_FPS) / dt))
    FRAME_STRIDE = max(1, FRAME_STRIDE)

    next_progress_t = 0.0

    # -----------------------------
    # Plot setup (°F)
    # -----------------------------
    fig, ax = plt.subplots()
    x_axis_in = np.linspace(0.0, L, n) / INtoM
    (line,) = ax.plot(x_axis_in, np.vectorize(K_to_F)(phi), linewidth=2)

    ax.grid(True)
    ax.set_xlim(0.0, L / INtoM)
    ax.set_xlabel("Depth (in)")
    ax.set_ylim(K_to_F(280.0), K_to_F(3200.0))
    ax.set_xlabel("Depth (m)")
    ax.set_ylabel("Temperature (°F)")

    frames = []
    melt_flag = False
    edge_log = []  # [t_s, wall_surface_F]

    # -----------------------------
    # Time loop
    # -----------------------------
    num_steps = int(round(tf / dt))

    for step in range(num_steps + 1):
        t = step * dt

        if t + 1e-12 >= next_progress_t:
            wallF = K_to_F(float(phi[0]))
            pct = 100.0 * t / tf if tf > 0 else 100.0
            print(f"Progress: t = {t:0.2f} / {tf:0.2f} s ({pct:0.1f}%) | Wall = {wallF:0.1f} °F")
            edge_log.append([t, wallF])
            next_progress_t += PROGRESS_DT

        Tw = phi[0]
        term_Ma = 1.0 + (gam - 1.0) / 2.0 * (Ma ** 2)

        sigma = 1.0 / (((0.5 * (Tw / T_g) * term_Ma + 0.5) ** 0.68) * (term_Ma ** 0.12))
        hG = hG_base * sigma

        d = (ap0 + 2.0 * DD) * np.ones(n)
        d[0] = (ap0 / 2.0) + DD + hG
        d[-1] = (ap0 / 2.0) + DD + h_nat

        a = -DD * np.ones(n - 1)
        b = -DD * np.ones(n - 1)

        c = ap0 * phi
        c[0] = (ap0 / 2.0) * phi[0] + hG * Taw
        c[-1] = (ap0 / 2.0) * phi[-1] + h_nat * T_amb

        phi = thomas_solve(a, d, b, c)

        if step % FRAME_STRIDE == 0:
            line.set_ydata(np.vectorize(K_to_F)(phi))
            ax.set_title(
            f"{config_path.stem} | Time: {t:.2f} s | "
            f"Pc: {Pc_chamber_psia:.0f} psia | "
            f"Wall: {K_to_F(float(phi[0])):.1f} °F"
        )
            fig.canvas.draw()
            buf = np.asarray(fig.canvas.buffer_rgba())
            frames.append(buf[:, :, :3].copy())

            if (phi[0] > T_melt) and (not melt_flag):
                print(f"WARNING: Wall surface melting at t = {t:.3f} s (Wall = {K_to_F(float(phi[0])):.1f} °F)")
                melt_flag = True

    plt.close(fig)

    # -----------------------------
    # Save outputs
    # -----------------------------
    imageio.mimsave(OUT_GIF, frames, duration=1.0 / GIF_FPS)

    with open(OUT_CSV, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["t_s", "wall_surface_F"])
        w.writerows(edge_log)

    sim_dt_per_frame = FRAME_STRIDE * dt
    print("Simulation complete.")
    print(f"GIF saved to: {OUT_GIF}")
    print(f"CSV saved to: {OUT_CSV}")
    print(f"Frames: {len(frames)} | sim-time per frame: {sim_dt_per_frame:.4f} s | GIF_FPS: {GIF_FPS}")


def main():
    # -----------------------------
    # CLI
    # -----------------------------
    if len(sys.argv) < 2:
        raise SystemExit("Usage: python heat_transfer_1d.py <config.json | folder_of_jsons>")

    targets = iter_config_paths(sys.argv[1])

    for config_path in targets:
        print(f"\n=== Running: {config_path.name} ===")
        try:
            run_one_config(config_path)
        except Exception as e:
            print(f"[FAIL] {config_path.name}: {e}")


if __name__ == "__main__":
    main()
