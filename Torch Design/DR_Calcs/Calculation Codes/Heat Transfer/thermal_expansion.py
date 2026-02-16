
import sys
import json
import pandas as pd
from pathlib import Path
from tkinter import Tk, filedialog

import matplotlib.pyplot as plt

IN2M = 0.0254

def open_args():

    # ---- FILE 1: Spreadsheet ----
    if len(sys.argv) >= 2:
        file_path = sys.argv[1]
    else:
        Tk().withdraw()
        file_path = filedialog.askopenfilename(
            title="Select spreadsheet file",
            filetypes=[("CSV files", "*.csv"),
                       ("Excel files", "*.xlsx *.xls")]
        )
        if not file_path:
            raise SystemExit("No spreadsheet selected.")

    try:
        if file_path.endswith(".csv"):
            df = pd.read_csv(file_path)
        elif file_path.endswith((".xlsx", ".xls")):
            df = pd.read_excel(file_path)
        else:
            raise SystemExit("Unsupported file type. Use .csv, .xlsx, or .xls")

        print("Successfully loaded Temperatures.\n")
        print(df.head())

    except Exception as e:
        raise SystemExit(f"Error reading spreadsheet: {e}")


    # ---- FILE 2: JSON config ----
    if len(sys.argv) >= 3:
        config_path = Path(sys.argv[2]).expanduser().resolve()
    else:
        Tk().withdraw()
        config_path = filedialog.askopenfilename(
            title="Select properties JSON file",
            filetypes=[("JSON files", "*.json")]
        )
        if not config_path:
            raise SystemExit("No JSON file selected.")

        config_path = Path(config_path).expanduser().resolve()

    try:
        with open(config_path, "r") as f:
            cfg = json.load(f)
    except Exception as e:
        raise SystemExit(f"Error reading JSON file: {e}")

    return df, cfg


def main():
    # check spreadsheet was given as argument and open it
    df, cfg = open_args()
    
    #unpack material properties and sim setup
    material_cfg = cfg["material"]
    areas_cfg = cfg["expansion_surfaces"]

    df["thermal_strain"] =  material_cfg["CTE"] * (df["wall_surface_F"] - material_cfg["CTE_ref_temp_F"])
    df["thermal_stress"] = df["thermal_strain"] * material_cfg["E_tension_Pa"]
    df["upper_force"] = df["thermal_stress"] * areas_cfg["A1_in"] * IN2M**2
    df["outer_force"] = df["thermal_stress"] * areas_cfg["A2_in"] * IN2M**2
    df["lower_force"] = df["thermal_stress"] * areas_cfg["A3_in"] * IN2M**2
    print(material_cfg["E_tension_Pa"])

    plt.plot(df["t_s"], df["wall_surface_F"])
    plt.show()
    print(df)

    return

if __name__ == "__main__":
     main()
