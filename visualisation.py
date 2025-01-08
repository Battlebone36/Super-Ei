import matplotlib.pyplot as plt
from present import Protein

class Visualise:
    def show_points(self, data: dict[tuple[int, int], tuple[str, int]]) ->  None:
        """Makes a plot of the protein."""
        # Filter the points out of the data
        # Polar
        x_p = [item[0][0] for item in data.items() if item[1][0] == "P"]
        y_p = [item[0][1] for item in data.items() if item[1][0] == "P"]

        # Hydrofobe
        x_h = [item[0][0] for item in data.items() if item[1][0] == "H"]
        y_h = [item[0][1] for item in data.items() if item[1][0] == "H"]

        # Cysteine        
        x_c = [item[0][0] for item in data.items() if item[1][0] == "C"]
        y_c = [item[0][1] for item in data.items() if item[1][0] == "C"]

        # Covalent bonds
        x_l = [item[0][0] for item in data.items()]
        y_l = [item[0][1] for item in data.items()]

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
