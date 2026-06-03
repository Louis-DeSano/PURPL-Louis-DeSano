# TCA Igniter Sizing Code
# 11/28/2025 : Created - Louis DeSano

# units standard: every variable is stored as SI, convert when inputting to external libraries if needed
# rocketCEA: uses English Engineering
# pyfluids: uses SI with Celsius

import numpy as np
import matplotlib.pyplot as plt
from rocketcea.cea_obj import CEA_Obj
from pyfluids import Fluid, FluidsList, Input
import pandas as pd

## Unit Conversions & Constants ##
n2lbf = 4.44822     # [N/lbf]

psi2Pa = 6894.76    # [psi/Pa]

lbm2kg = 0.453592   #[lbm/kg]

ft2m = 0.3048       # [ft/m]
m2in = 39.3701      # [m/in]

R2K = 0.555556      # [R/K]

# Constants
G0 = 9.81           # [m/s^2]
P_ATM = 101325      # [Pa]
T_AMB_CELSIUS = 20 # [deg C]

# Function: Calculate Area from Diameter
def diameter_to_area(diameter):
# Input: 
    # Diameter
# Output: 
    # Area
    return np.pi * (diameter/2)**2

# Function: Calculate Diameter from Area
def area_to_diameter(area):
# Input: 
    # Area
# Output: 
    # Diameter
    return 2 * np.sqrt(area / np.pi)

# Function: Calculate back pressure for choked flow
def choked_backpressure(k, stiffness, p_c):
# Inputs
    # k: specific heat ratio [--]
    # stiffness: amount stiffer than choked minimum as a decimal [--]
    # p_c: chamber pressure [Pa]

# Outputs
    # p_line: line pressure sized for choked flow [Pa]

    # calculate critical pressure ratio (minimum stiffness for choked flow)
    critRatio = (2 / (k+1)) ** (k / (k - 1)) 
    
    # add stiffness and calculate line pressure [Pa]
    p_line = (1 + stiffness) * critRatio**-1 * p_c

    return p_line

def unchoked_area(mDot, k, p_line, p_ratio, rho_line, Cd):
# Inputs
    # k: specific heat ratio [--]
    # stiffness: pressure drop / chamber pressure [--]
    # p_c: chamber pressure [Pa]

# Outputs
    # A_inj: orifice injection Area [m^2]

    A_inj = mDot/Cd * ( np.sqrt(2 * rho_line * p_line * (k / (k-1)) * (p_ratio**(2/k)  - p_ratio**((k+1) / k))) )**-1

    return A_inj

# Function: Construct line of excel sheet and add to sheet
def add_cad_param(sheet, name, dim, unit):
#Inputs
    #name: dimension name string
    #dim: dimension value
    #unit: unit string
#Output
    #sheet: Sheet with print to be printed to excel

    excelLine = [name, unit, f'{dim:0.3f}{unit}']
    sheet.append(excelLine)

    return sheet

# Function: Choked Flow at Flow Constriction
def throat_area(mDot, Cd, k, rho , p0):
# Inputs: 
    # k: Specific Heat Ratio    [--]
    # Cd: Discharge Coefficient [--]
    # rho: density              [kg/m^3]
    # p0: stagnation pressure   [Pa]
# Outputs:
    # At: Throat Area           [m^3]

    At = mDot / (Cd * np.sqrt(k * rho * p0 * ( 2 / (k + 1))**((k + 1)/(k - 1))))
    return At

# Function: Size Line for a velocity
def size_line(mDot, rho, v_line_chosen):
    # mDot: mass flow rate [kg/s]
    # rho: density [kg/m^3]
    # A_inj: injector area [m^2]
    # v_line_chosen: target line velocity [m/s]

    # Nominal SAE ORB outer diameters
    ID = np.arange(1/16, 2, 1/16) / m2in   # [m]

    # Line velocity for each ID:
    v_line = (mDot / rho) * (4 / (np.pi * ID**2))

    # Find the index of the diameter that gives closest velocity
    idx = np.abs(v_line - v_line_chosen).argmin()

    # Extract chosen diameter and velocity
    best_ID = ID[idx] * m2in
    best_v = v_line[idx]

    return best_ID, best_v

def main():

    ### Design Setpoints ###

    mDot_torch = 0.0408            # [kg/s] FROM MIE Study
    p_c = 300 * psi2Pa              # torch chamber pressure [psi->Pa]

    OF = 2                          # OF ratio

    Cd = 0.80                       # Discharge Coefficient
    choked = True

    # Propellant Specific
    ox_CEA = 'GOX'
    fuel_CEA = 'GH2'
    k_f = 1.41  # gh2 specific heat ratio
    k_ox = 1.40 # gox specific heat ratio  

    ### End Design Setpoints ###

    """ Mass Flow Rates """
    # calculate ox and fuel mass flow rates
    mDot_f = mDot_torch / (OF + 1) # [kg/s]
    mDot_ox = mDot_torch - mDot_f  # [kg/s]

    print('\nMass Flows')
    print(f' Torch Mass Flow [kg/s]: {mDot_torch:.5f}')
    print(f' Oxidizer Mass Flow [kg/s]: {mDot_ox:0.5f}')
    print(f' Fuel Mass Flow [kg/s]: {mDot_f:0.5f}')
    print('\n')


    """ Orifice/Throat Areas """
    # Calculate Line Pressures based on chosen stiffness (uniform stiffness)
    stiffness_ox = 0.35                    # deltaP / p_c [--]
    stiffness_f = 0.35                    # deltaP / p_c [--]


    #define pressure ratio and line pressure
    if (choked):
        p_line_ox = choked_backpressure(k_ox, stiffness_ox, p_c)
        p_line_f = choked_backpressure(k_f, stiffness_f, p_c)
    
    else:
        p_ratio_f = 1 / (1 + stiffness_f)   # pc/pline [--]
        p_line_f = p_c / p_ratio_f        # line pressure [Pa]
        p_ratio_ox = 1 / (1 + stiffness_ox)   # pc/pline [--]
        p_line_ox = p_c / p_ratio_ox          # line pressure [Pa]
    
    ## Get fluid properties to size orifice/throat areas ##
    #pyfluids for ox and fuel
    fuel_PYF = Fluid(FluidsList.Hydrogen).with_state(Input.pressure(p_line_f), Input.temperature(T_AMB_CELSIUS))
    ox_PYF = Fluid(FluidsList.Oxygen).with_state(Input.pressure(p_line_ox), Input.temperature(T_AMB_CELSIUS))

    rho_f = fuel_PYF.density            # [kg/m^3]
    rho_ox = ox_PYF.density             # [kg/m^3]
    sonic_f = fuel_PYF.sound_speed      # [m/s]
    sonic_ox = ox_PYF.sound_speed       # [m/s]

    print("\nPressures")
    print(f" Ox   [psi]: {p_line_ox/psi2Pa:0.3f}")
    print(f"   Stiffness: {stiffness_ox * 100}%")

    print(f" Fuel [psi]: {p_line_f/psi2Pa:0.3f}")
    print(f"   Stiffness: {stiffness_f * 100}%")

    print(f" Chamber [psi]: {p_c/psi2Pa:0.3f}")
    print("\n")

    # CEA for chamber exit (frozen gives lower performance bound)
    ispObj = CEA_Obj(oxName=ox_CEA, fuelName=fuel_CEA)
    isp = 0.8 * ispObj.estimate_Ambient_Isp(Pc= (p_c/psi2Pa), MR= OF, Pamb= (P_ATM/psi2Pa), eps=1, frozen=1)[0] # [s]
    rho_c = ispObj.get_Densities(Pc=(p_c/psi2Pa), MR= OF, eps=1, frozen=1)[0] * lbm2kg / (ft2m**3)             #[lbm/ft^3 -> kg/m^3]
    k_c = ispObj.get_Chamber_MolWt_gamma(Pc= (p_c/psi2Pa), MR= OF, eps=1)[1]

    thrust = mDot_torch * isp * G0                           # [N]
    cstar = ispObj.get_Cstar(Pc=p_c, MR=OF) * ft2m           # [m/s]
    T_comb = ispObj.get_Tcomb(Pc= p_c, MR= OF) * R2K         # [K]
    ## End of Fluid Properties ##

    # Assume Cd and cstar efficiency are unity to provide margin for choked flow
    A_t = throat_area(mDot_torch, Cd, k_c, rho_c, p_c)                  # Throat Area [m^2]
    #A_t = mDot_torch * cstar / (p_c)                                   # Throat Area [m^2]

    if choked:
        A_f = throat_area(mDot_f, Cd, k_f, rho_f, p_line_f)
        A_ox = throat_area(mDot_ox, Cd, k_ox, rho_ox, p_line_ox)
    else:
        A_f = unchoked_area(mDot_f, k_f, p_line_f, p_ratio_f, rho_f, Cd)        # Fuel Injection Area [m^2]
        A_ox = unchoked_area(mDot_ox, k_ox, p_line_ox, p_ratio_ox, rho_ox, Cd)    # Ox Injection Area   [m^2]

    D_f = area_to_diameter(A_f) 
    D_ox = area_to_diameter(A_ox)
    D_t = area_to_diameter(A_t)
    
    print("\nOrifice/Nozzle Sizing [in]")
    print(f" Fuel Injection Diameter: {D_f * m2in:0.3f}")
    print(f" Ox Injection Diameter: {D_ox* m2in:0.3f}")
    print(f" Throat Diameter: {D_t* m2in:0.3f}")
    print("\n")


    # Additional Outputs
    print("\nMisc")
    print(f" Thrust [N]: {thrust:0.3f} ")
    print(f" Tcomb  [K]: {T_comb:0.3f} ")
    print("\n")


    """Chamber Dimensions""" 
    #get chamber volume by defining stay time
    t_stay = 0.0005 # [s]
    V_chamber = t_stay * mDot_torch / rho_c # [m^3]
    Lstar = V_chamber / A_t # [m]
    #Lstar = 0.8 #m
    V_chamber = Lstar * A_t
    t_stay = V_chamber * rho_c / (mDot_torch)

    #chamber volume -> dimensions
    conv_angle = 65 # convergent angle [deg]
    r_contraction = 6.5 # contraction ratio

    A1 = r_contraction * A_t                                                # chamber area [m^2]
    D_c = area_to_diameter(A1)               
    L_conv = (D_c/2 - D_t/2) / np.sin(np.deg2rad(conv_angle))               # convergent length[m]
    L1 = ( V_chamber - A1*L_conv * (1 + np.sqrt(A_t/A1) + A_t/A1) ) / A1    # chamber length [m]

    print(f'\nChamber Sizing')
    print(f' Stay Time [s]: {t_stay:0.4f}')
    print(f' L Star [m]: {Lstar:0.3f}')
    print(f" Chamber Volume [m^3]: {V_chamber:0.5f}")
    print(f' Contraction Ratio: {A1/A_t:0.2f}')
    print(f' Convergent Angle [deg] {conv_angle}')
    print(f' Chamber Diameter [mm]: {D_c*10**3:0.3f}')
    print(f' Chamber Length [mm]: {L1*10**3:0.3f}')
    print(f' Convergent Length [mm]: {L_conv*10**3:0.3f}')
    print(f' Throat Diameter [mm]: {D_t*10**3:0.3f}')
    print("\n")


    """Line Sizing"""

    if choked:
        rho_star_ox = rho_ox * (2 / (k_ox+1)) ** (1 / (k_ox-1))
        rho_star_f = rho_f * (2 / (k_f+1)) ** (1 / (k_f-1))
        v_inj_ox = mDot_ox/(rho_star_ox*A_ox)
        v_inj_f = mDot_f/(rho_star_f*A_f)
    
    else:
        v_inj_ox = mDot_ox/(rho_ox*A_ox)
        v_inj_f = mDot_f/(rho_f*A_f)

    # find injeciton velocity
    D_line_ox, v_line_ox = size_line(mDot_ox, rho_ox, 50)
    D_line_f, v_line_f = size_line(mDot_f, rho_f, 160)

    """print("\nLine Sizing")
    print(f" Ox Line: ORB -{D_line_ox * 16:0.0f}")
    print(f" Ox Line Velocity [m/s]: {v_line_ox:0.3f}")
    print(f" Fuel Line: ORB -{D_line_f * 16:0.0f}")
    print(f" Fuel Line Velocity [m/s]: {v_line_f:0.3f}")"""

    print("\nOrifice/Nozzle Sizing [in]")
    print(f" Fuel Injection Diameter: {D_f * m2in:0.3f}")
    print(f" Ox Injection Diameter: {D_ox* m2in:0.3f}")
    print(f" Throat Diameter: {D_t* m2in:0.3f}")
    print("\n")

    print(f" Ox Injection Velocity [m/s]: {v_inj_ox:0.3f}")
    print(f" Fuel Injection Velocity [m/s]: {v_inj_f:0.3f}")
    print("\n")
    print(f"Ox Sonic [m/s]: {sonic_ox:0.3f}")
    print(f"Fuel Sonic [m/s]: {sonic_f:0.3f}")

    print(f'\nChamber Sizing IMPERIAL')
    print(f' Stay Time [s]: {t_stay:0.4f}')
    print(f' L Star [in]: {Lstar*m2in:0.3f}')
    print(f" Chamber Volume [in^3]: {V_chamber*m2in**3:0.5f}")
    print(f' Contraction Ratio: {A1/A_t:0.2f}')
    print(f' Convergent Angle [deg] {conv_angle}')
    print(f' Chamber Diameter [in]: {D_c*m2in:0.3f}')
    print(f' Chamber Length [in]: {L1*m2in:0.3f}')
    print(f' Convergent Length [in]: {L_conv*m2in:0.3f}')
    print(f' Throat Diameter [in]: {D_t*m2in:0.3f}')
    print("\n")

    print("Orifice/Nozzle Sizing [in]")
    print(f" Fuel Injection Diameter: {D_f * m2in:0.3f}")
    print(f" Ox Injection Diameter: {D_ox* m2in:0.3f}")
    print(f" Throat Diameter: {D_t* m2in:0.3f}")
    print("\n")

    print(f"fuel momentum {mDot_f*v_inj_f}")
    print(f"ox momentum {mDot_ox*v_inj_ox}")
    


    """Make Cad Parameter Sheet"""
    sheet = []

    sheet = add_cad_param(sheet, "Conv_Angle", conv_angle, "deg")

    sheet = add_cad_param(sheet, "d_t", D_t*10**3, "mm")
    sheet = add_cad_param(sheet, "d_ox", D_ox*10**3, "mm")
    sheet = add_cad_param(sheet, "d_f", D_f*10**3, "mm")
    sheet = add_cad_param(sheet, "d_c", D_c*10**3, "mm")

    sheet = add_cad_param(sheet, "L_c", L1*10**3, "mm")
    sheet = add_cad_param(sheet, "L_conv", L_conv*10**3, "mm")
    sheet = add_cad_param(sheet, "t_wall", 0.9, 'in')

    df = pd.DataFrame(sheet)
    df.to_csv("TCA_Igniter_dims.csv", index=False, header=False)
    
    return

if __name__== "__main__":
    main()
