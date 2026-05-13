import numpy as np
import itertools as it
import matplotlib.pyplot as plt

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
    "mission1": [4, 5, 4, 1, 2],
    "mission2": [3, 2, 3, 1, 1],
    "mission3": [1, 1, 3, 5, 5],
    "mission4": [2, 2, 2, 4, 3],
    "mission5": [2, 4, 3, 5, 4]
}
mission_scores = np.array(list(mission_scores_input.values()))  # type specified so no matrix multiplication issues later

def calculate_score(criteria, mission_scores):
    scores = mission_scores @ criteria
    winner = np.argmax(scores) + 1
    return scores, winner
# print(calculate_score(criteria, mission_scores))


delta_weight = np.arange(-0.1, 0.11, 0.05)
print(delta_weight)



def bar_chart_analysis(criteria, mission_scores):
    # getting all possibel weight combinations
    max_weight = 100 - (len(criteria) - 1) * 10
    combinations = it.product(range(10, max_weight + 1, 5), repeat=len(criteria))  # gets all permutations choosing 5 values from variation
    possible_weights = []
    for weights in combinations:
        if sum(weights) == 100:
            possible_weights.append(weights)

    possible_scores = np.empty((len(possible_weights), len(criteria)))
    possible_winners = []
    for i, weight in enumerate(possible_weights): 
        weight_array = np.array(weight).reshape(-1, 1) / 100

        scores, winner = calculate_score(weight_array, mission_scores) 
        possible_scores[i] = scores.ravel()
        possible_winners.append(winner)

    winner_count = np.bincount(possible_winners, minlength=len(mission_scores_input) + 1)[1:]

    scores_distribution = []
    for i in range(len(criteria)):
        scores_distribution.append(
            sum(possible_scores[:, i])
        )

    bins = range(1, len(criteria_input) + 1)
    # plt.bar(bins, scores_distribution)
    # plt.show()

    fig, ax1 = plt.subplots()

    # ax1.bar(bins, scores_distribution, alpha=0.5)
    ax1.set_xlabel("Mission")
    # ax1.set_ylabel("Cummulative score across all weight combinations")

    ax2 = ax1.twinx()
    ax2.bar(bins, winner_count)
    ax2.set_ylabel("Number of wins")

    plt.xticks(bins)
    plt.title("Sensitivity Analysis")
    plt.show()
    return scores_distribution, winner_count

# bar_chart_analysis(criteria, mission_scores)