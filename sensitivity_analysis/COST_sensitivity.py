'''COST SENSITIVITY ANALYSIS
    1) edit the equations used to determine the cost -- LINE 13
    2) edit what the uncertainty you want to analyse is, this is a +/- percentage -- LINE 38
    3) edit the base costs without the singular equations etc -- LINE 44
    4) edit the base cost without all the equations -- LINE 53
    5) CHECK #VARIABLES -- LINES 55, 56, 66
    6) replace the names for the chart -- LINE 66'''

from itertools import product
import numpy as np
import matplotlib.pyplot as plt

# equations used to determine the cost 
def cost_eq1(uncertainty):
    cost = 0.1*1 + 100
    return cost * uncertainty

def cost_eq2(uncertainty):
    cost = 0.2*1**2 + 200
    return cost * uncertainty 

def cost_eq3(uncertainty):
    cost = 0.3*1**3 + 30
    return cost * uncertainty 

cost_equations = [cost_eq1, cost_eq2, cost_eq3]

# now finding what that equation output would be through full uncertainty range
def vary_cost_eq(cost_eq, relative_uncertainty):
    cost_variation = []
    uncertainty_variation = np.linspace(-relative_uncertainty, relative_uncertainty, 10*relative_uncertainty)
    for uncertainty in uncertainty_variation:
        cost = cost_eq(uncertainty/100)
        cost_variation.append(cost)
    return cost_variation

cost_vars = []
relative_uncertainty = 10  # PARAMETER TO CHANGE
for cost_eq in cost_equations:
    cost_variation = vary_cost_eq(cost_eq, relative_uncertainty)
    cost_vars.append(cost_variation)

cost_sums = []
base_costs = [180, 130, 190]  # base cost without ONE of the equations
# variation of ONE equation
for i, base_cost in enumerate(base_costs):
    cost_sum = []
    for var in cost_vars[i]:
        cost_sum.append(base_cost + var)
    cost_sums.append(cost_sum)

# getting the variation of EVERY equation
base_cost = 110  # base cost without ANY of the equations
cost_sum = []
for var1, var2, var3 in product(cost_vars[0], cost_vars[1], cost_vars[2]):
    cost_sum.append(base_cost + var1 + var2 + var3)
cost_sums.append(cost_sum)

for cost_sum in cost_sums:
    min_cost = min(cost_sum)
    max_cost = max(cost_sum)
    avg_cost = sum(cost_sum)/len(cost_sum)
    print("===========================")
    print(min_cost, avg_cost, max_cost)

costs = {'all': cost_sums[-1], 'eq1': cost_sums[0], 'eq2': cost_sums[1], 'eq3': cost_sums[2]}
plt.boxplot(costs.values(), labels=costs.keys(), orientation='horizontal')
# plt.boxplot(cost_sums, orientation='horizontal')
plt.xlabel(f"Cost [million €]")
plt.show()