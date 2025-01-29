import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from code.classes.protein import Protein
import numpy as np
from timeit import default_timer as timer
import pandas as pd
import seaborn as sns


# class Visualise:
def filter_data(data: dict[tuple[int, int, int], tuple[str, int]], amino: str) -> tuple[list[int], list[int], list[int]]:
    """
    Helper function for visualise_protein().
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
    Visualizes the 3D structure of a protein using matplotlib using filter_data().
    This function generates a 3D plot of the protein structure, where different types of amino acids 
    (Polar, Hydrophobic, Cysteine) are represented by different colors. It also marks non-sequential 
    bonds with different colors and line styles based on their bond type.
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
    Tests the given algorithm 100 times and collects performance data.

    This function runs the provided algorithm on a predefined protein sequence
    100 times, recording the stability scores and execution times for each run.
    It then calculates and returns the maximum stability score, minimum stability
    score.
    """
    sequence = "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP"
    test_protein = Protein(sequence)
    stability_scores: list[int] = []
    time_scores: list[float] = []

    # Test the algorithm 100 times and store the result
    for i in range(100):
        print(i)
        start = timer()
        alg = algorithm(test_protein)
        result_protein = alg.run()
        end = timer()
        stability_scores.append(result_protein.stability())
        time_scores.append(end - start)

    # Define variables for making the plot
    max_score = max(stability_scores)
    min_score = min(stability_scores)
    return (stability_scores, max_score, min_score)


def plot_algorithm_split(algorithm, ax, width: float, offset: float) -> None:
    """
    Helper function
    Plot the gathered data as a histogram in the total plot with bars split.

    Uses:
    - hist_of_algorithm()
    """
    # Define variables
    stability_scores, max_score, min_score = hist_of_algorithm(algorithm)

    # Define locations and heights of the bars
    bins = np.arange(min_score, max_score + 1)
    hist_values, bin_edges = np.histogram(stability_scores, bins=bins)
    x_positions = (bin_edges[:-1] + bin_edges[1:]) / 2 + offset + width

    # Make plot
    ax.bar(x_positions,
           hist_values,
           width=width,
           edgecolor="black",
           label=f"{algorithm.__name__}")


def plot_algorithm_together(algorithm, ax, color: str) -> None:
    """
    Helper function
    Plot the gathered data as a histogram in the total plot with bars together.
    
    Uses:
    - hist_of_algorithm()
    """
    # Define variables
    stability_scores, max_score, min_score = hist_of_algorithm(algorithm)

    # Make plot with bars together
    ax.hist(
        stability_scores,
        bins= max_score - min_score + 1,
        range=(min_score - 0.5, max_score + 0.5),
        density=True,
        label=f"{algorithm.__name__}",
        linewidth=5,
        histtype="step"
    )

def visualise_algorithm(algorithms, command: str="split") -> None:
    """
    Visualise the distribution of stability scores of the given algorithms. 
    The plot style can be altered using the command parameter.

    command:
    ----------------
    - "split" (default): gives a histogram with the bars split.
    - "together": gives a histogram with the bars together.
    - "line": gives a histogram in a line plot.

    Uses:
    ----------------
    - plot_algorithm_split()
    - plot_algorithm_together()
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
    else:
        raise ValueError

    # Plot settings
    plt.xlabel("Stability")
    plt.ylabel("Occurrency")
    plt.title("Functionality of algorithm out of 100 tests")
    plt.legend()
    plt.show()


def import_algorithm_data(path, df):
    """
    Import the data out of a certain algorithm csv files.
    """
    df_data_gathered = pd.read_csv(path)
    return pd.concat([df, df_data_gathered], axis=0)


def make_plot(df, algorithms, axes=None, type="many_runs") -> None:
    """
    Generate various types of plots to visualize the performance of different algorithms.

    Parameters:
    -----------
    df (pd.DataFrame): The dataframe containing the data to be plotted. It must include columns 'Algorithm', 'Stability', 'Time', and 'Iteration'.
    algorithms (list): A list of algorithm functions whose performance is to be visualized.
    axes (list, optional): A list of matplotlib axes objects for plotting subplots. Defaults to None.
    type (str, optional): The type of plot to generate. Options are 'many_runs', 'time', and 'iteration_stability'. Defaults to 'many_runs'.

    Plot Types:
    ------------
    - 'many_runs': Generates histograms showing the stability distribution of each algorithm.
    - 'time': Generates a boxplot showing the run time distribution of each algorithm.
    - 'iteration_stability': Generates line plots showing the stability over iterations for each algorithm.

    Example:
    --------
    make_plot(df, [algorithm1, algorithm2], axes=axes, type="many_runs")
    """
    
    # Make a histogram by default and by type occurency-stability
    if type == "many_runs":
        for i, algorithm in enumerate(algorithms):
            bins = df["Stability"].max() - df["Stability"].min()
            df_filtered = df[df["Algorithm"] == f"{algorithm.__name__}"]
            print(f"{algorithm.__name__:<20} mean: {df_filtered['Stability'].mean():<25} min: {df_filtered['Stability'].min()}")
            plot = sns.histplot(
                df_filtered["Stability"],
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
            plot.set_ylim(0, 0.16)
        plot.set_ylabel("Chance")

    # If type is time-stability make a line plot with mean time on y-axis
    # and stability on x-axis
    elif type == "time":
        order = []
        for algorithm in algorithms:
            order.append(f"{algorithm.__name__}")
            df_filtered = df[df["Algorithm"] == f"{algorithm.__name__}"]
            print(f"{algorithm.__name__:<20} mean: {df_filtered['Time'].mean():<20} min: {df_filtered['Time'].min():<20} max: {df_filtered['Time'].max()}")
        # print(df.head(20))
        plot = sns.boxplot(
            data=df,
            x="Algorithm",
            y="Time",
            order=order
        )
        plot.set_ylabel("Time (s)")
        plot.set_ylim(0, 30)
        plot.set_title("Run Time of the Algorithms")

    elif type == "iteration_stability":
        # Plot every algorithm        
        for i, algorithm in enumerate(algorithms):

            # Filter the data per algorithm
            if algorithm.__name__ == "HillClimb":
                df_filtered = df[(df["Algorithm"] == f"{algorithm.__name__}") & (df["Iteration"] <= 4500)] # Data after 4500 is not succesfull
                df_filtered_4500 = df_filtered[df_filtered["Iteration"] == 4500]
                mean = df_filtered_4500["Stability"].mean()
                minimum = df_filtered_4500["Stability"].min()
                print(f"{algorithm.__name__:<20} mean: {mean:<20} min: {minimum:<20}")
            else:
                df_filtered = df[df["Algorithm"] == f"{algorithm.__name__}"]
                df_filtered_4998 = df_filtered[df_filtered["Iteration"] == 4998]
                mean = df_filtered_4998["Stability"].mean()
                minimum = df_filtered_4998["Stability"].min()

                print(f"{algorithm.__name__:<20} mean: {mean:<20} min: {minimum:<20}")
            

            # Make the multiple plots or single plot
            if len(algorithms) > 1:
                plot = sns.lineplot(
                    data=df_filtered,
                    x="Iteration",
                    y="Stability",
                    ax=axes[i]
                )
            else:
                plot = sns.lineplot(
                    data=df_filtered,
                    x="Iteration",
                    y="Stability"
                )

            # Set the plot settings
            plot.set_title(f"{algorithm.__name__}")
            plot.set_xlabel("Iterations") 
            plot.set_ylabel("Stability")
            plot.set_xlim(0, 1000)
            plot.set_ylim(-22, 0)
            plot.set_yticks(np.arange(-22, 1, 2))

            # Verbose
            print(f"{algorithm.__name__} is done")
        plot.legend(["Mean Stability", "95% Confidence Interval"])
    plt.show()

def visualise_algorithm_data(algorithms, type: str="many_runs") -> None:
    """
    Imports data and makes plots for the given algorithms.

    Parameters:
    ----------------
    algorithms (list): A list of algorithm classes to visualize.
    type (str): The type of visualization to create. Options are:
        - "iteration_stability": Visualizes the stability of iterations.
        - "many_runs": Visualizes the results of many runs.
        - "time": Visualizes the time-based performance.
        - "one_hour_run": Visualizes the results of a one-hour run.
    """
    # Make DataFrame and import data from files
    type = type.lower()
    df = pd.DataFrame()
    if type == "iteration_stability":
        for algorithm in algorithms:
            path = f"code/data/iteration_stability/{algorithm.__name__}.csv"
            df = import_algorithm_data(path=path, df=df)
    elif type == "many_runs" or type == "time":
        for algorithm in algorithms:
            path = f"code/data/many_runs/{algorithm.__name__}.csv"
            df = import_algorithm_data(path=path, df=df)
    elif type == "one_hour_run":
        for algorithm in algorithms:
            path = f"code/data/one_hour_run/{algorithm.__name__}.csv"
            df = import_algorithm_data(path=path, df=df)

    
    # Make specific plot
    if type == "many_runs" or type == "one_hour_run":
        fig, axes = plt.subplots(1, len(algorithms), figsize=(15, 5), sharex=True, sharey=True)
        make_plot(df=df, algorithms=algorithms, axes=axes)
    elif type == "time":
        grouped = df.groupby(["Algorithm", "Protein", "Stability"]).mean().reset_index()
        make_plot(df=grouped, algorithms=algorithms, type=type)
    elif type == "iteration_stability":
        fig, axes = plt.subplots(1, len(algorithms), figsize=(15, 5), sharex=True, sharey=True)
        make_plot(df=df, algorithms=algorithms, axes=axes, type=type)


