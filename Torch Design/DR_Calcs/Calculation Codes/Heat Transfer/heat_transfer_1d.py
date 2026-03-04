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
def K_to_F(Tk: float) -> float: return (Tk - 273.15) * 9.0 / 5.0 + 32.0
def R_to_K(TR: float) -> float: return TR * (5.0 / 9.0)
def psi_to_Pa(psi: float) -> float: return psi * 6894.757293168

def get_script_dir() -> Path:
    try: return Path(__file__).resolve().parent
    except NameError: return Path.cwd()

def load_materials():
    mat_path = get_script_dir() / "materials.json"
    if not mat_path.exists():
        raise FileNotFoundError("Could not find materials.json in the script directory.")
    with open(mat_path, "r") as f:
        return json.load(f)

def load_config(config_path: Path) -> dict:
    with open(config_path, "r") as f:
        cfg = json.load(f)
    return cfg

def iter_config_paths(arg: str):
    p = Path(arg).expanduser().resolve()
    if p.is_file(): return [p], p.parent
    if p.is_dir():
        cfgs = sorted(p.glob("*.json"))
        # Filter out materials.json so we don't try to run it as a config
        cfgs = [c for c in cfgs if c.name != "materials.json"]
        if not cfgs: raise SystemExit(f"No configs found in: {p}")
        return cfgs, p
    raise SystemExit(f"Path not found: {p}")

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
# Simulation Logic
# -----------------------------
class SimulationResult:
    def __init__(self, name, x_in, t_melt_f, frames_data, times, melt_time, edge_log):
        self.name, self.x_in, self.t_melt_f = name, x_in, t_melt_f
        self.frames_data, self.times, self.melt_time, self.edge_log = frames_data, times, melt_time, edge_log

def run_simulation(config_path: Path, material_db: dict):
    cfg = load_config(config_path)
    cea_cfg, wall_cfg = cfg["cea"], cfg["wall"]
    geom_cfg, solv_cfg = cfg["geometry"], cfg["solver"]

    # Pull material properties from database
    mat_name = wall_cfg.get("material")
    if mat_name not in material_db:
        raise ValueError(f"Material '{mat_name}' not found in materials.json")
    
    mat = material_db[mat_name]
    k_w = float(mat["k_w"])
    rho = float(mat["rho"])
    cp_wall = float(mat["cp_wall"])
    T_melt = float(mat["T_melt_K"])
    h_nat = float(mat["h_nat"])
    T_amb = float(mat["T_amb_K"])

    # CEA & Physics
    Pc_psia = float(cea_cfg["Pc_chamber_psia"])
    MR, eps, frozen = float(cea_cfg["MR"]), float(cea_cfg.get("eps", 1.0)), int(cea_cfg.get("frozen", 0))
    D_in = float(geom_cfg["D_in"])
    A = np.pi * (D_in * INtoM / 2.0)**2
    Ma = float(geom_cfg.get("Ma", 1.0))
    Thr_D = float(geom_cfg.get("Throat_D_in", D_in)) * INtoM
    Thr_A = np.pi * (Thr_D / 2.0)**2
    
    L, n = float(wall_cfg["L_in"]) * INtoM, int(wall_cfg["n"])
    tf, dt = float(solv_cfg["tf"]), float(solv_cfg["dt"])
    stride = max(1, int(round((1.0 / int(cfg.get("output", {}).get("gif_fps", 30))) / dt)))

    cea = CEA_Obj(oxName=cea_cfg["oxName"], fuelName=cea_cfg["fuelName"])
    T_ch_R = cea.get_Temperatures(Pc=Pc_psia, MR=MR, eps=eps, frozen=frozen)[0]
    T_g = R_to_K(float(T_ch_R))
    _, gam = cea.get_Chamber_MolWt_gamma(Pc=Pc_psia, MR=MR, eps=eps)
    Cp_c, visc_c, cond_c, Pr = cea.get_Chamber_Transport(Pc=Pc_psia, MR=MR, eps=eps, frozen=frozen)
    cp_g, mu, k_g = float(Cp_c)*4184.0, float(visc_c)*1e-4, float(cond_c)*0.4184
    Pc_Pa, Cstar = psi_to_Pa(Pc_psia), float(cea.get_Cstar(Pc=Pc_psia, MR=MR)) * 0.3048

    # Pre-calc
    dy = L / n
    phi = float(wall_cfg.get("T_init_K", T_amb)) * np.ones(n)
    hG_base = (0.026/(Thr_D**0.2)) * ((mu**0.2 * cp_g)/(Pr**0.6)) * (Pc_Pa/Cstar)**0.8 * (Thr_A/A)**0.9
    Taw = T_g * (1.0 + (gam - 1.0)/2.0 * Ma**2 * Pr**(1/3))
    DD, ap0 = k_w/dy, rho*cp_wall*dy/dt
    
    history, times, edge_log = [], [], []
    melt_time = None

    # Time Loop
    for step in range(int(tf/dt) + 1):
        t = step * dt
        sigma = 1.0 / (((0.5*(phi[0]/T_g)*(1.0 + (gam-1.0)/2.0*Ma**2) + 0.5)**0.68) * ((1.0 + (gam-1.0)/2.0*Ma**2)**0.12))
        hG = hG_base * sigma

        d = (ap0 + 2.0*DD)*np.ones(n); d[0] = ap0/2 + DD + hG; d[-1] = ap0/2 + DD + h_nat
        a, b = -DD*np.ones(n-1), -DD*np.ones(n-1)
        c = ap0*phi; c[0] = (ap0/2)*phi[0] + hG*Taw; c[-1] = (ap0/2)*phi[-1] + h_nat*T_amb
        phi = thomas_solve(a, d, b, c)

        if melt_time is None and phi[0] >= T_melt: melt_time = t
        if step % stride == 0:
            history.append(phi.copy())
            times.append(t)
            edge_log.append([t, K_to_F(phi[0])])
            
    return SimulationResult(f"{config_path.stem} ({mat_name})", np.linspace(0, L, n)/INtoM, K_to_F(T_melt), history, times, melt_time, edge_log)

def save_individual_data(res, folder):
    folder.mkdir(parents=True, exist_ok=True)
    with open(folder / f"{res.name.split(' ')[0]}_log.csv", "w", newline="") as f:
        csv.writer(f).writerows([["t_s", "wall_surface_F"]] + res.edge_log)

def create_gif(results, output_path, fps=30):
    fig, ax = plt.subplots(figsize=(10, 6))
    lines = []
    ax.set_xlim(0, max(r.x_in[-1] for r in results))
    ax.set_ylim(K_to_F(280), 3500); ax.set_xlabel("Depth (in)"); ax.set_ylabel("Temperature (°F)"); ax.grid(True)

    for r in results:
        m_lbl = f" (Melt: {r.melt_time:.2f}s)" if r.melt_time else " (No Melt)"
        l, = ax.plot(r.x_in, K_to_F(r.frames_data[0]), label=f"{r.name}{m_lbl}", linewidth=2)
        ax.axhline(y=r.t_melt_f, color=l.get_color(), linestyle='--', alpha=0.4)
        lines.append(l)
    
    ax.legend(loc='upper right', fontsize='small')
    frames = []
    for i in range(max(len(r.frames_data) for r in results)):
        for idx, r in enumerate(results):
            lines[idx].set_ydata(K_to_F(r.frames_data[min(i, len(r.frames_data)-1)]))
        ax.set_title(f"Thermal Profile | Time: {results[0].times[min(i, len(results[0].times)-1)]:.2f} s")
        fig.canvas.draw()
        frames.append(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy())

    imageio.mimsave(output_path, frames, fps=fps)
    plt.close(fig)

def main():
    if len(sys.argv) < 2: raise SystemExit("Usage: python script.py <config_or_folder>")
    targets, input_root = iter_config_paths(sys.argv[1])
    material_db = load_materials()
    results_base = get_script_dir() / f"{input_root.name} Results"
    results_base.mkdir(parents=True, exist_ok=True)

    mode = "individual"
    if input_root.is_dir() and len(targets) > 1:
        if input("Run together? (t/i): ").lower() == 't': mode = "together"

    all_res = []
    for cfg_p in targets:
        print(f"Solving: {cfg_p.name}")
        res = run_simulation(cfg_p, material_db)
        indiv_folder = results_base / cfg_p.stem
        save_individual_data(res, indiv_folder)
        create_gif([res], indiv_folder / f"{cfg_p.stem}.gif")
        if mode == "together": all_res.append(res)

    if mode == "together" and all_res:
        create_gif(all_res, results_base / "Combined_Analysis.gif")

if __name__ == "__main__":
    main()