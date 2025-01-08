import matplotlib.pyplot as plt
import numpy as np


class Protein:
    
    def __init__(self, sequence: str) -> None:
        sequence = sequence.upper()
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
    
    def give_data(self) -> dict[tuple[int, int], tuple[str, int]]:
        return self.data

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
        borders = [min(x_l), max(x_l), min(y_l), max(y_l)]

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
        plt.xlim(min(borders) - 1, max(borders) + 1)
        plt.ylim(min(borders) - 1, max(borders) + 1)
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
        score = 0
        for acid in h_acids.items():
            friends = self.neighbours(acid[0])
            for friend in friends:
                
                if friend in self.data and abs(acid[1][1] - self.data[friend][1]) == 1:
                    continue
                elif self.h_bond(acid[0], friend):
                    score -= 1
        score //= 2
                    
        print(score)
        

        return 1
    
    # Doesn't work yet but makes code cleaner
    def rotate_coord(pivot: tuple[int, int], coord: tuple[int, int], matrix) -> tuple[int, int]:
        """Rotates the coordinate around a pivot with a matrix"""
        rel_coord = (coord[0] - pivot[0], coord[1] - pivot[1])
        v_rel_coord = np.array(rel_coord)
        rel_rot_coord = matrix @ v_rel_coord
        rel_new_coord = tuple(rel_rot_coord)
        new_coord = (rel_new_coord[0] + pivot[0], rel_new_coord[1] + pivot[1])
        return new_coord

    def fold(self, pivot: tuple[int, int], direction: str) -> bool:
        # Raise an error if the coordinate is not in the dataset
        if pivot not in self.data:
            print("coordinate not in dataset")
            raise(IndexError)
        
        # Make the command case insensitive and define the rotation matrix
        direction = direction.lower()
        if direction == "right":
            rot_matrix = np.array(self.right_turn)
        elif direction == "left":
            rot_matrix = np.array(self.left_turn)
        
        # Store the points that need to be rotated
        points_to_rot: dict[tuple[int, int], tuple[str, int]] = {}
        for acid in self.data.items():
            if acid[1][1] > self.data[pivot][1]:
                points_to_rot[acid[0]] = acid[1]
        
        # The function above rotate_coord doesn't work yet
        # Check if the rotation doesn't clash with the existing folding
        for acid in points_to_rot:
            rel_coord = (acid[0] - pivot[0], acid[1] - pivot[1])
            v_rel_coord = np.array(rel_coord)
            rel_rot_coord = rot_matrix @ v_rel_coord
            rel_new_coord = tuple(rel_rot_coord)
            new_coord = (rel_new_coord[0] + pivot[0], rel_new_coord[1] + pivot[1])
            if new_coord in self.data:
                return False
        
        # Rotate every point
        for acid in points_to_rot:
            rel_coord = (acid[0] - pivot[0], acid[1] - pivot[1])
            v_rel_coord = np.array(rel_coord)
            rel_rot_coord = rot_matrix @ v_rel_coord
            rel_new_coord = tuple(rel_rot_coord)
            new_coord = (rel_new_coord[0] + pivot[0], rel_new_coord[1] + pivot[1])
            old_point = self.data.pop(acid)
            self.data[new_coord] = old_point
        return True








protein1 = Protein("HHPHPC")
protein1.stability()
protein1.show_points()
protein1.fold((2, 0), "left")
protein1.fold((2, 2), "right")

