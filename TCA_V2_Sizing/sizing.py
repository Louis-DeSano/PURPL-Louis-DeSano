import pint as pt
import numpy as np
import cea 

import ceaPrinter

ureg = pt.UnitRegistry()

Pc_range = np.array([15,60]) * ureg.bar
thrust_range = np.array([1000, 6000]) * ureg.lbf
OF_range = np.array([1,5])

reac_names = ["C3H8O,2propanol", "O2(L)"]
T_reactant = np.array([298, 90]) * ureg.K
fuel_weights = np.array([1.0, 0.0])
oxidant_weights = np.array([0.0, 1.0])

reac = cea.Mixture(reac_names)
prod = cea.Mixture(reac_names, products_from_reactants=True)

solver = cea.RocketSolver(prod, reactants=reac)
solution = cea.RocketSolution(solver)

OF = 2.1
Pc = 500 * ureg.psi
Patm = 1 * ureg.bar
Pratio = Pc/Patm
supar = [25, 50, 75]

weights = reac.of_ratio_to_weights(oxidant_weights, fuel_weights, OF)
hc = reac.calc_property(cea.ENTHALPY, weights, T_reactant.to(ureg.K).magnitude)/cea.R

solver.solve(solution, weights, Pc.to(ureg.bar).magnitude, Pratio.magnitude, supar=supar, hc=hc, iac=True)
ceaPrinter.print_full_ouput(solution)





