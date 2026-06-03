from rocketcea.cea_obj import CEA_Obj, add_new_fuel
import numpy as np
import matplotlib.pyplot as plt

# n-Butanol (1-butanol) card
# h,cal = heat of formation = -159,600 cal/mol (liquid, 298.15 K)
# rho = 0.810 g/cc
card_str = """
fuel  C4H10O(L)  C 4  H 10  O 1
      wt%=100.00  h,cal=-78196.0  t(k)=298.15  rho=0.810
"""

add_new_fuel('nButanol', card_str)

R2F = -458.67 # rankine to farenheit

OF_sweep = np.arange(1,3,0.01)


fuel_list = ["Methanol", "Ethanol", "nButanol", "Isopropanol", "JetA"]
ox = "LOX"

pc = 500 #psi

for fuel in fuel_list:
    temp_array = []
    cstar_array = []
    isp_array = []
    for OF in OF_sweep:
        
        ispObj = CEA_Obj(fuelName=fuel, oxName=ox)
        temp_array.append(ispObj.get_Tcomb(pc,OF) + R2F)
        isp_array.append(ispObj.get_Isp(pc,OF))
        cstar_array.append(ispObj.get_Cstar(pc,OF))
    
    plt.figure(1)
    plt.subplot(2,1,1)
    plt.plot(OF_sweep,temp_array)
    plt.grid(True)
    plt.legend(fuel_list)
    plt.ylabel("Combustion Temperature [F]")
    plt.title("Combustion Temperature vs. OF Ratio")

    plt.subplot(2,1,2)
    plt.plot(OF_sweep,cstar_array)
    plt.grid(True)
    plt.legend(fuel_list)
    plt.title("c* vs. OF Ratio")
    plt.xlabel("OF Ratio")
    plt.ylabel("c* [ft/s]")

    plt.figure(2)
    plt.plot(temp_array, cstar_array)
    plt.grid(True)
    plt.legend(fuel_list)
    plt.xlabel("Combustion Temperature [F]")
    plt.ylabel("Cstar [ft/s]")
    plt.title("Cstar vs. Temperature")



plt.figure(1)
plt.subplot(2,1,1)
ax = plt.gca()
ax_r = ax.twinx()
y1, y2 = ax.get_ylim()
ax_r.set_ylim((y1-32)*5/9 + 273.15, (y2-32)*5/9 + 273.15)
ax_r.set_ylabel("Combustion Temperature [K]")
plt.subplot(2,1,2)
ax = plt.gca()
ax_r = ax.twinx()
y1, y2 = ax.get_ylim()
ax_r.set_ylim(y1*0.3048, y2*0.3048)
ax_r.set_ylabel("c* [m/s]")

plt.figure(2)
# right y-axis (ft/s → m/s)
ax = plt.gca()
ax_r = ax.twinx()
y1, y2 = ax.get_ylim()
ax_r.set_ylim(y1*0.3048, y2*0.3048)
ax_r.set_ylabel("c* [m/s]")
# top x-axis (°F → K)
ax_t = ax.twiny()
x1, x2 = ax.get_xlim()
ax_t.set_xlim((x1-32)*5/9 + 273.15, (x2-32)*5/9 + 273.15)
ax_t.set_xlabel("Combustion Temperature [K]")

plt.show()