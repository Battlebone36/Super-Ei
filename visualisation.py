import matplotlib.pyplot as plt
from protein import Protein

class Visualise:
    def filter_data(data: dict[tuple[int, int], tuple[str, int]], amino: str = None) -> tuple[list[int], list[int]]:
        if amino is None:
            x = [item[0][0] for item in data.items()]
            y = [item[0][1] for item in data.items()]
        else:
            x = [item[0][0] for item in data.items() if item[1][0] == amino]
            y = [item[0][1] for item in data.items() if item[1][0] == amino]
        return (x, y)

    def visualise_protein(protein: Protein) ->  None:
        """Makes a plot of the protein."""
        # Filter the points out of the data
        data: dict[tuple[int, int], tuple[str, int]] = protein.give_data()
        x_p, y_p = Visualise.filter_data(data, "P")
        x_h, y_h = Visualise.filter_data(data, "H")
        x_c, y_c = Visualise.filter_data(data, "C")
        x_l, y_l = Visualise.filter_data(data)

        # Borders for plot and legend
        borders = [min(x_l), max(x_l), min(y_l), max(y_l)]

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
        fig, ax = plt.subplots()
        ax.plot(x_l, y_l, c = "black", alpha= 0.8, linewidth= 5)
        ax.plot(x_p, y_p, "bo", markersize= 20)
        ax.plot(x_h, y_h, "ro", markersize= 20)
        ax.plot(x_c, y_c, "go", markersize= 20)
        plt.xlim(min(borders) - 1, max(borders) + 1)
        plt.ylim(min(borders) - 1, max(borders) + 1)
        plt.legend(in_plot)
        plt.show()

protein_vis = Protein("H")
Visualise.visualise_protein(protein_vis)
