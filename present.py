import matplotlib.pyplot as plt

class Protein:
    
    def __init__(self, sequence: str) -> None:
        self.data: dict[tuple[int, int], tuple[str, int]] = {}
        for i, char in enumerate(sequence):
            self.data[(i, 0)] = (f"{char}", i)

    def show_points(self) ->  None:
        # Filter the points out of the data
        # Polair
        x_p = [item[0][0] for item in self.data.items() if item[1][0] == "P"]
        y_p = [item[0][1] for item in self.data.items() if item[1][0] == "P"]

        # Hydrofobe
        x_h = [item[0][0] for item in self.data.items() if item[1][0] == "H"]
        y_h = [item[0][1] for item in self.data.items() if item[1][0] == "H"]

        # Cystine        
        x_c = [item[0][0] for item in self.data.items() if item[1][0] == "C"]
        y_c = [item[0][1] for item in self.data.items() if item[1][0] == "C"]

        # Covalent bonds
        x_l = [item[0][0] for item in self.data.items()]
        y_l = [item[0][1] for item in self.data.items()]

        # Borders for plot and legend
        x_min, x_max, y_min, y_max = min(x_l), max(x_l), min(y_l), max(y_l)
        in_plot = []
        if x_l:
            in_plot.append("Bond")
        if x_p:
            in_plot.append("Polair")
        if x_h:
            in_plot.append("Hydrofobe")
        if x_c:
            in_plot.append("Cystine")

        # Make the plot with dots and line
        fig, ax = plt.subplots()
        ax.plot(x_l, y_l, c = "black", alpha= 0.8, linewidth= 5)
        ax.plot(x_p, y_p, "bo", markersize= 20)
        ax.plot(x_h, y_h, "ro", markersize= 20)
        ax.plot(x_c, y_c, "go", markersize= 20)
        plt.xlim(x_min - 1, x_max + 1)
        plt.ylim(y_min - 1, y_max + 1)
        plt.legend(in_plot)
        plt.show()

    def options(self, x: int, y: int) -> None:
        "Prints the the North, East, South and West coördinates"
        x_diff = [0, 1, 0, -1]
        y_diff = [1, 0, -1, 0]
        result = [(x_diff[i] + x, y_diff[i] + y) for i in range(4)]
        print(result)


protein1 = Protein("HHPHPC")
# protein1.options(5, 5)
protein1.show_points()


# def calculate_stability(data: list[list[str | tuple[str, int]]]) -> int:
#     filtered_list = list(filter(lambda x: x[0] == "H", filtered_list))
#     filtered_dict = {}
#     score = 0
#     for value in filtered_list:
#         filtered_dict[f"{value[2], value[3]}"] = value[0]
#     for value in filtered_list:
#         for option in options(value[2], value[3]):
#             if f"{option}" in filtered_dict and filtered_dict[f"{option}"] == "H":
#                 print("ja")
#     return 1
