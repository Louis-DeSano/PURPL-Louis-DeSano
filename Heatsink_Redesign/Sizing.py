from rocketcea.cea_obj import CEA_Obj
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
PSI2PA = 6894.76
FT2M = 0.3048
LBM2KG = 0.453592
N2LBF = 0.224809


# Setpoints (Imperial)
mdot_array = [20, 22, 24]   # lbm/s

Pc = 500.0 #psi
Pe = 12.0 #psi
Pamb = 14.7 #psi

eta_cstar = 0.875
eta_cF = 0.875

fuel = "Isopropanol"
ox = "LOX"

OF_sweep = np.arange(1.0, 3.5, 0.01)

# SI Conversion

Pc_si = Pc * PSI2PA
Pe_si = Pe * PSI2PA
Pamb_si = Pamb * PSI2PA

Pc_psi = Pc
Pe_psi = Pe
Pamb_psi = Pamb

cea = CEA_Obj(fuelName=fuel, oxName=ox)

plt.figure(figsize=(8,5))

for mdot in mdot_array:

    mdot_si = mdot * LBM2KG

    F_list = []

    for OF in OF_sweep:

        cstar_ft = cea.get_Cstar(Pc=Pc_psi, MR=OF)

        eps = cea.get_eps_at_PcOvPe(Pc=Pc_psi,MR=OF,PcOvPe=Pc_psi/Pe_psi)

        cf = cea.get_PambCf(Pc=Pc_psi, MR=OF,Pamb=Pamb_psi,eps=eps)[0]
        cf *= eta_cF
        cstar_si = cstar_ft * FT2M * eta_cstar

        At_si = mdot_si * cstar_si / Pc_si

        F_si = cf * Pc_si * At_si

        F_lbf = F_si * N2LBF

        F_list.append(F_lbf)

    F_array = np.array(F_list)

    plt.plot(OF_sweep, F_array, label=f"ṁ = {mdot} lbm/s")

plt.xlabel("O/F Ratio")
plt.ylabel("Thrust (lbf)")
plt.title("LOX / Isopropanol Thrust vs Mixture Ratio\nPc=500 psi, Pe=12 psi")

plt.grid()

plt.gca().yaxis.set_major_formatter(StrMethodFormatter('{x:,.0f}'))

plt.legend()

plt.show()