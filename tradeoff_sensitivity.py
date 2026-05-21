import numpy as np
import itertools as it
import matplotlib.pyplot as plt

'''HOW TO USE
    1) adjust the criteria_input with your different criteria and their current weights
    2) adjust the mission_scores_input with the score each mission receives for each criteria
    RUN THE CODE!
    the sensitivity analysis basically varies the weights for all possible values with 5% increments with min weight of 10, and max weight dont worry about it (it's 
    max weight possible where every other criteria has a weight of 10)
    it then shows you the number of times each mission won for each weight variation
    then graphs it'''

# adjust this as needed
criteria_input = {
    "lca": 0.3,
    "collision": 0.45,
    "prop": 0.25
}
criteria = np.array(list(criteria_input.values())).reshape(-1, 1) 

# here you put in the scores for each category (in order of the criteria)
# if it's in excel, just copy paste this and your table into chat and ask it to write it out
mission_scores_input = {
    "STR": [3, 4, 2],
    "MTR": [4, 3, 2],
    "PSW": [2, 2, 2],
    "M&T": [2, 3, 4],
    "ISC": [2, 3, 5]
}
mission_scores = np.array(list(mission_scores_input.values()))  # type specified so no matrix multiplication issues later

def calculate_score(criteria, mission_scores):
    scores = mission_scores @ criteria
    max_score = np.max(scores)
    winner = np.flatnonzero(np.isclose(scores, max_score))  # accounting for a tie situation (assigning a win to all who tied)
    return winner
# print(calculate_score(criteria, mission_scores))

# delta_weight = np.arange(-0.1, 0.11, 0.05)
# print(delta_weight)

def bar_chart_analysis(criteria, mission_scores):
    # getting all possibel weight combinations
    max_weight = 100 - (len(criteria) - 1) * 10
    combinations = it.product(it.chain([0], range(10, max_weight + 1, 5)), repeat=len(criteria))  # gets all permutations choosing 5 values from variation
    possible_weights = []
    for weights in combinations:
        if sum(weights) == 100:
            possible_weights.append(weights)

    winner_count = np.zeros(len(mission_scores_input))

    for weight in possible_weights: 
        weight_array = np.array(weight).reshape(-1, 1) / 100
        winners = calculate_score(weight_array, mission_scores) 
        
        for winner in winners:
            winner_count[winner] += 1
    
    plt.bar(['STR', 'MTR', 'PSW', 'M&T', 'ISC'], winner_count/len(possible_weights))
    plt.ylabel("Proportion of wins")
    plt.xlabel("Mission concept")
    plt.title("Sustainability sensitivity analysis")
    plt.ylim(0, 1.1)
    plt.show()
    return winner_count

bar_chart_analysis(criteria, mission_scores)


'''CAPTURE MECHANISM'''
# criteria_input = {
#     "mass": 0.2,
#     "tumbling_performance": 0.1,
#     "docking_performance": 0.25,
#     "operational_flexibility": 0.1,
#     "reliability_and_redundancy": 0.25,
#     "TRL": 0.1
# }
# mission_scores_input = {
#     "RM":       [3, 4, 5, 5, 5, 4],
#     "GAG":      [5, 3, 4, 4, 3, 4],
#     "Clamp":    [4, 3, 4, 4, 5, 4],
#     "PA":       [4, 3, 4, 3, 5, 5],
#     "MEV-M":    [5, 3, 4, 3, 5, 5]
# }

'''CHEMICAL PROPULSION'''
# criteria_input = {
#     "isp": 0.2,
#     "thurst": 0.1,
#     "storage": 0.25,
#     "TRL": 0.1
# }
# mission_scores_input = {
#     "LMP-103S":     [3, 5, 5, 5],
#     "ASCENT/ADN":   [2, 3, 4, 5],
#     "SHP163/HAN":   [4, 2, 4, 5],
#     "HTP-kerosene": [4, 5, 4, 3]
# }

'''ELECTRIC PROPULSION'''
# criteria_input = {
#     "isp": 0.15,
#     "thurst": 0.3,
#     "dry_mass": 0.2,
#     "power": 0.2,
#     "TRL": 0.15
# }
# mission_scores_input = {
#     "NEXT-C":       [3, 5, 2, 2, 5],
#     "PPS-1350E":    [3, 4, 4, 2, 5],
#     "PPS-1350S":    [3, 3, 4, 3, 5],
#     "PPS-5000":     [3, 5, 3, 2, 5],
#     "PPS-X00":      [3, 3, 5, 3, 4],
#     "Solar Sail":   [5, 1, 3, 5, 3]
# }

'''SUSTAINABILITY'''
# criteria_input = {
#     "lca": 0.35,
#     "collision_probability": 0.25,
#     "propulsion_efficiency": 0.4
# }
# mission_scores_input = {
#     "STR":  [3, 4, 2],
#     "MTR":  [4, 2, 2],
#     "PSW":  [2, 4, 2],
#     "M&T":  [2, 3, 4],
#     "ISC":  [2, 2, 5]
# }

'''FEASIBILITY & RELIABILITY'''
# criteria_input = {
#     "cost": 0.25,
#     "complexity": 0.15,
#     "schedule": 0.2,
#     "redundancy": 0.15,
#     "risk": 0.25
# }
# mission_scores_input = {
#     "STR":  [4, 5, 4, 1, 2],
#     "MTR":  [3, 2, 3, 1, 1],
#     "PSW":  [1, 1, 3, 5, 5],
#     "M&T":  [2, 2, 3, 4, 3],
#     "ISC":  [2, 4, 3, 5, 4]
# }

'''INNOVATION'''
# criteria_input = {
#     "team_preference": 0.4,
#     "technical_originality": 0.25,
#     "future_potential": 0.35
# }
# mission_scores_input = {
#     "STR":  [1, 2, 3],
#     "MTR":  [1, 3, 4],
#     "PSW":  [1, 4, 4],
#     "M&T":  [5, 4, 4],
#     "ISC":  [1, 2, 3]
# }

'''TECHNICAL PERFORMANCE'''
# criteria_input = {
#     "time": 0.4,
#     "equivalent_mass": 0.6
# }
# mission_scores_input = {
#     "STR":  [3, 3],
#     "MTR":  [5, 2],
#     "PSW":  [3, 4],
#     "M&T":  [3, 5],
#     "ISC":  [4, 1]
# }

'''MISSION CONCEPT'''
# criteria_input = {
#     "innovation": 0.15,
#     "technical_performance": 0.35,
#     "sustainability": 0.2,
#     "feasibility_and_reliability": 0.3
# }
# mission_scores_input = {
#     "STR":  [1.95, 3.6, 2.85, 3.2],
#     "MTR":  [2.55, 3.2, 2.7, 2.05],
#     "PSW":  [2.8, 3.6, 2.5, 3],
#     "M&T":  [4.4, 4.2, 3.05, 2.75],
#     "ISC":  [1.95, 2.2, 3.2, 3.45]
# }