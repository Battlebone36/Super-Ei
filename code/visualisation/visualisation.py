import matplotlib.pyplot as plt
from code.classes.protein import Protein


# class Visualise:
def filter_data(data: dict[tuple[int, int], tuple[str, int]], amino: str | None = None) -> tuple[list[int], list[int], list[int]]:
    """"Filters the coordinates of amino acids from a dataset, optionally filtering by type."""
    if amino is None:
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
    x_l, y_l, z_l = filter_data(data)

    # # Borders for plot and legend
    # borders = [min(x_l), max(x_l), min(y_l), max(y_l), min(z_l), max(z_l)]

    in_plot = []
    if x_l:
        in_plot.append("Bond")
    if x_p:
        in_plot.append("Polar")
    if x_h:
        in_plot.append("Hydrofobe")
    if x_c:
        in_plot.append("Cysteine")

    # Make the plot with dots and line
    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    ax.plot(x_l, y_l, z_l, c = "black", alpha=0.8, linewidth=4)
    ax.plot(x_p, y_p, z_p, "bo", markersize=12)
    ax.plot(x_h, y_h, z_h, "ro", markersize=12)
    ax.plot(x_c, y_c, z_c, "go", markersize=12)

    # Mark non-sequential bonds
    line_info = {-1: ("r", 2), -5: ("g", 3.5)}
    for acid in data.items():
        friends = protein.neighbours(acid[0])
        for friend in friends:
            if friend in data and abs(acid[1][1] - data[friend][1]) != 1:
                score = protein.type_bond(acid[0], friend)
                if score < 0:
                    ax.plot([acid[0][0], friend[0]], [acid[0][1], friend[1]], line_info[score][0], linewidth=line_info[score][1], linestyle="dotted")

    # Specific plot settings
    ax.set_aspect("equal", adjustable="box")
    ax.grid(False)
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])

    ax.set_axis_off()

    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    ax.xaxis.pane.set_edgecolor((1.0, 1.0, 1.0, 0.0))
    ax.yaxis.pane.set_edgecolor((1.0, 1.0, 1.0, 0.0))
    ax.zaxis.pane.set_edgecolor((1.0, 1.0, 1.0, 0.0))

    # Show plot
    stability = protein.stability()
    plt.title(f"Stability score: {stability}", fontsize=17, fontweight="bold")
    plt.xlim(min(x_l) - 1, max(x_l) + 1)
    plt.ylim(min(y_l) - 1, max(y_l) + 1)
    plt.legend(in_plot)
    plt.show()
