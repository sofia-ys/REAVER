import numpy as np
import itertools as it

# adjust this as needed
criteria_input = {
    "cost": 0.25,
    "complexity": 0.20,
    "schedule": 0.20,
    "redundancy": 0.10,
    "risk": 0.25
}
criteria = np.array(list(criteria_input.values())).reshape(-1, 1) 

# here you put in the scores for each category (in order of the criteria)
# if it's in excel, just copy paste this and your table into chat and ask it to write it out
mission_scores_input = {
    "mission1": [4, 5, 3, 2, 4],
    "mission2": [1, 2, 3, 1, 2],
    "mission3": [3, 1, 1, 5, 2],
    "mission4": [2, 2, 3, 4, 3],
    "mission5": [1, 3, 3, 4, 3]
}
mission_scores = np.array(list(mission_scores_input.values()))  # type specified so no matrix multiplication issues later

def calculate_score(criteria, mission_scores):
    scores = mission_scores @ criteria
    winner = np.argmax(scores)
    return scores.flatten(), winner
print(calculate_score(criteria, mission_scores))

# defining the ranges of possible weights
    
# def sensitivity_analysis(criteria):
#     return 0

# getting all possibel weight combinations
max_weight = 100 - (len(criteria) - 1) * 10
combinations = it.product(range(10, max_weight + 1, 5), repeat=5)  # gets all permutations choosing 5 values from variation
possible_weights = []
for weights in combinations:
    if sum(weights) == 100:
        possible_weights.append(weights)

possible_scores = []
for weight in possible_weights:
    scores, _ = calculate_score(weight, mission_scores)
    possible_scores.append(scores)

    print("-------------------------------------------------------------------\n")
    print(weight)
    print(scores)
    # print(mission_scores)
# print(possible_scores)
