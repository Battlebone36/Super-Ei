import matplotlib.pyplot as plt

class Protein:
    
    def __init__(self, sequence: str) -> None:
        self.data: dict[tuple[int, int], tuple[str, int]] = {}
        # for i, char in enumerate(sequence):
        #     self.data[(i, 0)] = (f"{char}", i)

        # test data for a random protein
        self.data[(0, 0)]  = ("P", 0)
        self.data[(1, 0)]  = ("H", 1)
        self.data[(1, 1)]  = ("H", 2)
        self.data[(0, 1)]  = ("H", 3)
        self.data[(0, 2)]  = ("H", 4)
        self.data[(1, 2)]  = ("P", 5)
        self.data[(2, 2)]  = ("P", 6)
        self.data[(2, 1)]  = ("H", 0)
        self.data[(2, 0)]  = ("P", 0)
        self.data[(2, -1)] = ("H", 0)
        self.data[(1, -1)] = ("P", 0)
        self.data[(0, -1)] = ("P", 0)


        self.left_turn = [[0, -1], [1, 0]]
        self.right_turn = [[0, 1], [-1, 0]]

    def show_points(self) ->  None:
        """Makes a plot of the protein"""
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

    def neighbours(self, coord: tuple[int, int]) -> list[tuple[int, int]]:
        "Prints the North, East, South and West coördinates of the one given"
        x_diff = [0, 1, 0, -1]
        y_diff = [1, 0, -1, 0]
        result = [(x_diff[i] + coord[0], y_diff[i] + coord[1]) for i in range(4)]
        return result
    
    def h_bond(self, coord1: tuple[int, int], coord2: tuple[int, int]) -> bool:
        if coord1 in self.data and coord2 in self.data:
            return self.data[coord1][0] == "H" and self.data[coord2][0] == "H"
        return False


    def stability(self) -> int:
        # Filter out the "H" acids
        h_acids: dict[tuple[int, int], tuple[str, int]] = {}
        for acid in self.data.items():
            if acid[1][0] == "H":
                h_acids[acid[0]] = acid[1]
        print(h_acids)
                
        #Loop throug the "H" acids and look to the neighbours
        for acid in h_acids.items():
            friends = self.neighbours(acid[0])
            for friend in friends:
                print(self.h_bond(acid[0], friend))

        score = 0

        return 1



protein1 = Protein("HHPHPC")
protein1.stability()
protein1.show_points()

