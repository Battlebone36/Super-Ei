import matplotlib.pyplot as plt
from code.classes.protein import Protein
from code.algorithms.randomise import random_fold
import numpy as np
from timeit import default_timer as timer
import pandas as pd
import seaborn as sns


# class Visualise:
def filter_data(data: dict[tuple[int, int, int], tuple[str, int]], amino: str) -> tuple[list[int], list[int], list[int]]:
    """"
    Filters the coordinates of amino acids from a dataset, optionally filtering by type.
    """
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
    """
    Makes a plot of the protein.
    """
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
    labeled_scores = set()
    for acid in data.items():
        friends = protein.neighbours(acid[0])
        for friend in friends:
            if friend in data and abs(acid[1][1] - data[friend][1]) != 1:
                score = protein.type_bond(acid[0], friend)
                if score in line_info:
                    # Create a label for in the legend
                    label = f"{score} bond" if score not in labeled_scores else ""
                    ax.plot([acid[0][0], friend[0]], 
                            [acid[0][1], friend[1]], 
                            [acid[0][2], friend[2]], 
                            line_info[score][0], 
                            linewidth=line_info[score][1], 
                            linestyle="dotted",
                            label=label)
                    labeled_scores.add(score)

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
    for i in range(1000):
        print(i)
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
    Helper function
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
    ax.bar(x_positions,
           hist_values,
           width=width,
           edgecolor="black",
           label=f"{algorithm.__name__}, mean: {mean_time:0.2}s ± {std_time:0.2}s")
    # , mean: {mean_time:0.2}s ± {std_time:0.2}s


def plot_algorithm_together(algorithm, ax, color: str) -> None:
    """
    Helper function
    Plot the gathered data as a histogram in the total plot with bars together.
    
    Uses:
    - hist_of_algorithm()
    """
    # Define variables
    stability_scores, max_score, min_score, mean_time, std_time = hist_of_algorithm(algorithm)

    # Make plot with bars together
    ax.hist(
        stability_scores,
        # edgecolor="black",
        bins= max_score - min_score + 1,
        range=(min_score - 0.5, max_score + 0.5),
        density=True,
        # alpha=0.8,
        label=f"{algorithm.__name__}",
        linewidth=5,
        histtype="step"
    )
    # , mean: {mean_time:0.2}s ± {std_time:0.2}s
    # histtype="barstacked",
    
            # density=True


def plot_algorithm_line(algorithm, ax) -> None:
    """
    Helper function
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
    ax.plot(x_positions,
            hist_values,
            label=f"{algorithm.__name__}")
    # , mean: {mean_time:0.2}s ± {std_time:0.2}s


def visualise_algorithm(algorithms, command: str="split") -> None:
    """
    Visualise the distribution of stability scores of the given algorithms. 
    The plot style can be altered using the command parameter.

    command:
    - "split" (default): gives a histogram with the bars split.
    - "together": gives a histogram with the bars together.
    - "line": gives a histogram in a line plot.

    Uses:
    - plot_algorithm_split()
    - plot_algorithm_together()
    - plot_algorithm_line()
    """
    # Define total plot and dataframe
    fig = plt.figure()
    ax = fig.add_subplot(111)
    colors = {
        "random_fold": "blue",
        "greedy_search_sequence": "green",
        "climbing_fold": "red",
        "beter_climbing_fold": "orange",
        "even_beter_climbing_fold": "yellow" 
    }

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


def import_algorithm_data(algorithm, df):
    """
    Import the data out of a certain algorithm csv files.
    """

    df_data_gathered = pd.read_csv(f"code/data/csv_data/{algorithm.__name__}.csv")
    return pd.concat([df, df_data_gathered], axis=0)

def make_plot(df, algorithms, axes=None, type="occurency-stability") -> None:
    """
    Make the a histogram of the occurency-stability of all the algorithms in the dataset
    """
    
    # Make a histogram by default and by type occurency-stability
    if type == "occurency-stability":
        for i, algorithm in enumerate(algorithms):
            bins = df["stability"].max() - df["stability"].min()
            df_filtered = df[df["algorithm"] == f"{algorithm.__name__}"]
            plot = sns.histplot(
                df_filtered["stability"],
                kde=True,
                bins=bins,
                stat="density",
                discrete=True,
                ax=axes[i]
            )
            ticks = np.arange(0, -bins - 1, -bins // 8)
            plot.set_xticks(ticks=ticks)
            plot.set_title(f"{algorithm.__name__}")
            plot.text(-bins, 0.95, s=f'Trials: {len(df_filtered)}', fontsize=10, color='black')
        plot.set_xlabel("Stability")
        plot.set_ylabel("Chance")
        plot.set_ylim(0, 1)
        

    # If type is time-stability make a line plot with mean time on y-axis
    # and stability on x-axis
    elif type == "time-stability":
        plot = sns.lineplot(
            data=df,
            x="stability",
            y="time",
            hue="algorithm",
            marker="o",
        )
        plot.set_xlabel("Stability")
        plot.set_ylabel("Mean time")
        plot.set_xlim(right=1)
        plot.legend()

    plt.show()


def visualise_algorithm_data(algorithms, type: str="occurency-stability") -> None:
    """
    Import all data and make the plot for the data.
    """
    # Make DataFrame and import data from files
    type = type.lower()
    df = pd.DataFrame()
    for algorithm in algorithms:
        df = import_algorithm_data(algorithm=algorithm, df=df)
    
    # Make specific plot
    if type == "occurency-stability":
        fig, axes = plt.subplots(1, len(algorithms), figsize=(15, 5), sharex=True, sharey=True)
        make_plot(df=df, algorithms=algorithms, axes=axes)
    elif type == "time-stability":
        grouped = df.groupby(["algorithm", "protein", "stability"]).mean().reset_index()
        make_plot(df=grouped, algorithms=algorithms, type=type)
