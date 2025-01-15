import matplotlib.pyplot as plt
from code.classes.protein import Protein
from code.algorithms.randomise import random_fold
import numpy as np
from timeit import default_timer as timer


# class Visualise:
def filter_data(data: dict[tuple[int, int, int], tuple[str, int]], amino: str) -> tuple[list[int], list[int], list[int]]:
    """"Filters the coordinates of amino acids from a dataset, optionally filtering by type."""
    if amino == "all":
        x = [item[0][0] for item in data.items()]
        y = [item[0][1] for item in data.items()]
        z = [item[0][2] for item in data.items()]
    else:
        x = [item[0][0] for item in data.items() if item[1][0] == amino]
        y = [item[0][1] for item in data.items() if item[1][0] == amino]
        z = [item[0][2] for item in data.items() if item[1][0] == amino]
    return (x, y, z)

def visualise_protein(protein: Protein) -> None:
    """Makes a plot of the protein."""
    # Filter the points out of the data
    data: dict[tuple[int, int, int], tuple[str, int]] = protein.give_data()
    x_p, y_p, z_p = filter_data(data, "P")
    x_h, y_h, z_h = filter_data(data, "H")
    x_c, y_c, z_c = filter_data(data, "C")
    x_l, y_l, z_l = filter_data(data, "all")

    # Make the plot with dots and line
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(x_l, y_l, z_l, "black", markersize=12, linewidth=4, label="Bond")
    ax.plot(x_p, y_p, z_p, "bo", markersize=12, linewidth=4, label="Polar")
    ax.plot(x_h, y_h, z_h, "ro", markersize=12, linewidth=4, label="Hydrofobe")
    ax.plot(x_c, y_c, z_c, "go", markersize=12, linewidth=4, label="Cysteine")

    # Mark non-sequential bonds
    line_info = {-1: ("r", 2), -5: ("g", 3.5)}
    for acid in data.items():
        friends = protein.neighbours(acid[0])
        for friend in friends:
            if friend in data and abs(acid[1][1] - data[friend][1]) != 1:
                score = protein.type_bond(acid[0], friend)
                if score < 0:
                    ax.plot([acid[0][0], friend[0]], [acid[0][1], friend[1]], [acid[0][2], friend[2]], line_info[score][0], linewidth=line_info[score][1], linestyle="dotted")

    # Make the distance between the points equal and remove the axis
    ax.set_aspect("equal", adjustable="box")
    ax.set_axis_off()

    # Show plot
    stability = protein.stability()
    plt.title(f"Stability score: {stability}", fontsize=17, fontweight="bold")
    plt.legend()
    plt.show()

def hist_of_algorithm(algorithm) -> tuple[list[int], int, int, float, float]:
    """
    2nd generation helper function

    Test the algorithm 100 times and store the data into lists.
    Return the data with some calculated time variables.
    """
    sequence = "HHPHHHPHPHHHPH"
    test_protein = Protein(sequence)
    stability_scores: list[int] = []
    time_scores: list[float] = []

    # Test the algorithm 100 times and store the result
    for i in range(100):
        start = timer()
        result_protein: Protein = algorithm(test_protein)
        end = timer()
        stability_scores.append(result_protein.stability())
        time_scores.append(end - start)

    # Define variables for making the plot
    max_score = max(stability_scores)
    min_score = min(stability_scores)
    mean_time = np.mean(time_scores)
    std_time = np.std(time_scores)
    return (stability_scores, max_score, min_score, mean_time, std_time)


def plot_algorithm_split(algorithm, ax, width: float, offset: float) -> None:
    """
    Plot the gathered data as a histogram in the total plot with bars split.

    Uses:
    - hist_of_algorithm()
    """
    # Define variables
    stability_scores, max_score, min_score, mean_time, std_time = hist_of_algorithm(algorithm)

    # Define locations and heights of the bars
    bins = np.arange(min_score, max_score + 1)
    hist_values, bin_edges = np.histogram(stability_scores, bins=bins)
    x_positions = (bin_edges[:-1] + bin_edges[1:]) / 2 + offset + width

    # Make plot
    ax.bar(x_positions, hist_values, width=width, edgecolor="black",
           label=f"{algorithm.__name__}, mean: {mean_time:0.2}s ± {std_time:0.2}s")


def plot_algorithm_together(algorithm, ax) -> None:
    """
    Helper function.
    Plot the gathered data as a histogram in the total plot with bars together.
    
    Uses:
    - hist_of_algorithm()
    """
    # Define variables
    stability_scores, max_score, min_score, mean_time, std_time = hist_of_algorithm(algorithm)

    # Make plot with bars together
    ax.hist(stability_scores, edgecolor="black", bins= max_score - min_score + 1,
            range=(min_score - 0.5, max_score + 0.5), alpha=0.8, label=f"{algorithm.__name__}, mean: {mean_time:0.2}s ± {std_time:0.2}s")


def plot_algorithm_line(algorithm, ax) -> None:
    """
    Helper function.
    Plot the gathered data as histogram in a line plot.

    Uses:
    - hist_of_algorithm()
    """
    # Define variables
    stability_scores, max_score, min_score, mean_time, std_time = hist_of_algorithm(algorithm)
    
    # Define locations and heights of the line
    bins = np.arange(min_score, max_score + 1)
    hist_values, bin_edges = np.histogram(stability_scores, bins=bins)
    x_positions = (bin_edges[:-1] + bin_edges[1:]) / 2 + 0.5

    # Make the line plot
    ax.plot(x_positions, hist_values, label=f"{algorithm.__name__}, mean: {mean_time:0.2}s ± {std_time:0.2}s")


def visualise_algorithm(algorithms, command: str="split") -> None:
    """
    Visualise all algorithms given. The plot can be altered with a command.

    command:
    - "split" (default) gives a histogram with the bars split
    - "together" gives a histogram with the bars together
    - "line"     gives a histogram in a line plot.

    Uses:
    - plot_algorithm_split()
    - plot_algorithm_together()
    - plot_algorithm_line()
    """
    # Define total plot and dataframe
    fig = plt.figure()
    ax = fig.add_subplot(111)

    # Check commands
    if command == "split":
        bar_width = 1 / (len(algorithms) + 1)
        for i, algorithm in enumerate(algorithms):
            plot_algorithm_split(algorithm=algorithm, ax=ax, width=bar_width, offset=bar_width * i)
    elif command == "together":
        for algorithm in algorithms:
            plot_algorithm_together(algorithm=algorithm, ax=ax)
    elif command == "line":
        for algorithm in algorithms:
            plot_algorithm_line(algorithm=algorithm, ax=ax)
    else:
        raise ValueError

    # Plot settings
    plt.xlabel("Stability")
    plt.ylabel("Occurrency")
    plt.title("Functionality of algorithm out of 100 tests")
    plt.legend()
    plt.show()