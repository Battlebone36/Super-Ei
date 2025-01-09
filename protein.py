import numpy as np


class Protein:
    
    def __init__(self, sequence: str) -> None:
        """
        Initialise the protein from a sequence of "P", "H" and "C" characters.
        The character with its order are stored in a dictionary where the keys are the coordinates.
        """
        sequence = sequence.upper()
        self.data: dict[tuple[int, int], tuple[str, int]] = {}
        # for i, char in enumerate(sequence):
        #     self.data[(i, 0)] = (f"{char}", i)

        # test data for a random protein
        self.data[(0, 0)]  = ("P", 0)
        self.data[(1, 0)]  = ("C", 1)
        self.data[(1, 1)]  = ("H", 2)
        self.data[(0, 1)]  = ("H", 3)
        self.data[(0, 2)]  = ("H", 4)
        self.data[(1, 2)]  = ("P", 5)
        self.data[(2, 2)]  = ("P", 6)
        self.data[(2, 1)]  = ("H", 7)
        self.data[(2, 0)]  = ("C", 8)
        self.data[(2, -1)] = ("H", 9)
        self.data[(1, -1)] = ("C", 10)
        self.data[(0, -1)] = ("P", 11)


        self.left_turn = [[0, -1], [1, 0]]
        self.right_turn = [[0, 1], [-1, 0]]
    
    def give_data(self) -> dict[tuple[int, int], tuple[str, int]]:
        return self.data

    def neighbours(self, coord: tuple[int, int]) -> list[tuple[int, int]]:
        """Returns the North, East, South and West coördinates of the one given."""
        x_diff = [0, 1, 0, -1]
        y_diff = [1, 0, -1, 0]
        result = [(x_diff[i] + coord[0], y_diff[i] + coord[1]) for i in range(4)]
        return result
    
    def h_bond(self, coord1: tuple[int, int], coord2: tuple[int, int]) -> bool:
        """
        Returns if the coordinates are in the protein dataset and
        if they have a bond that makes the protein stronger.
        """
        if coord1 in self.data and coord2 in self.data:
            return self.data[coord1][0] == "H" and self.data[coord2][0] == "H"
        return False

    def hc_bond(self, coord1: tuple[int, int], coord2: tuple[int, int]) -> bool:
        if coord1 in self.data and coord2 in self.data:
            return ((self.data[coord1][0] == "C" and self.data[coord2][0] == "H") or
                    (self.data[coord1][0] == "H" and self.data[coord2][0] == "C"))
        return False
    
    def c_bond(self, coord1: tuple[int, int], coord2: tuple[int, int]) -> bool:
        if coord1 in self.data and coord2 in self.data:
            return self.data[coord1][0] == "C" and self.data[coord2][0] == "C"
        return False

    def stability(self) -> int:
        """A function that calculates the stability of a protein and returns it in an integer."""
        # Filter out the "H" acids
        h_acids: dict[tuple[int, int], tuple[str, int]] = {}
        for acid in self.data.items():
            if acid[1][0] == "H" or acid[1][0] == "C":
                h_acids[acid[0]] = acid[1]
                
        # Loop through the "H" acids and look to the neighbours
        score = 0
        for acid in h_acids.items():
            friends = self.neighbours(acid[0])
            for friend in friends:
                
                if friend in self.data and abs(acid[1][1] - self.data[friend][1]) == 1:
                    continue
                elif self.h_bond(acid[0], friend):
                    score -= 1
                elif self.hc_bond(acid[0], friend):
                    score -= 1
                elif self.c_bond(acid[0], friend):
                    score -= 5
        score //= 2
        return score
    
    # Doesn't work yet but makes code cleaner
    def rotate_coord(self, coord: tuple[int, int], pivot: tuple[int, int],  matrix) -> tuple[int, int]:
        """
        Rotates the coordinate around a pivot with a rotation matrix.

        Input:
        - coord: the point that is being rotated.
        - pivot: is the point around which the coord is rotated.
        - matrix: the rotation matrix that defines the direction of the rotation.

        Output:
        - Rotated coordinate        
        """
        rel_coord = (coord[0] - pivot[0], coord[1] - pivot[1])
        v_rel_coord = np.array(rel_coord)
        rel_rot_coord = matrix @ v_rel_coord
        rel_new_coord = tuple(rel_rot_coord)
        new_coord = (rel_new_coord[0] + pivot[0], rel_new_coord[1] + pivot[1])
        return new_coord

    def fold(self, pivot: tuple[int, int], direction: str) -> bool:
        """Folds a protein at a given pivot point in a certain direction: "left" or "right"."""
        # Raise an error if the coordinate is not in the dataset.
        if pivot not in self.data:
            print("coordinate not in dataset")
            raise(IndexError)
        
        # Make the command case insensitive and define the rotation matrix.
        direction = direction.lower()
        if direction == "right":
            rot_matrix = np.array(self.right_turn)
        elif direction == "left":
            rot_matrix = np.array(self.left_turn)
        
        # Store the points that are following in the sequence.
        # These are following in the sequence
        points_to_rot: dict[tuple[int, int], tuple[str, int]] = {}
        for acid in self.data.items():
            if acid[1][1] > self.data[pivot][1]:
                points_to_rot[acid[0]] = acid[1]
        
        # Check if the rotation doesn't clash with the existing folding.
        for acid in points_to_rot:
            new_coord = self.rotate_coord(acid, pivot, rot_matrix)
            if new_coord in self.data:
                return False
        
        # Rotate every point
        for acid in points_to_rot:
            new_coord = self.rotate_coord(acid, pivot, rot_matrix)
            old_point = self.data.pop(acid)
            self.data[new_coord] = old_point
        return True



protein1 = Protein("HHPHPC")
