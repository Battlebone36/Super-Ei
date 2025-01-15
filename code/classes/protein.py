import numpy as np

class Protein:
    def __init__(self, sequence: str, command: str | None = None) -> None:
        """
        Initialise the protein from a sequence of "P", "H" and "C" characters.
        The character with its sequence index are stored in a dictionary where the keys are the coordinates.
        """
        self.sequence = sequence.upper()
        self.data: dict[tuple[int, int, int], tuple[str, int]] = {}
        self.load_data(command)
        self.rotations = {
            "x_pos": np.array([[1, 0, 0], [0, 0, 1], [0, -1, 0]]),
            "x_neg": np.array([[1, 0, 0], [0, 0, -1], [0, 1, 0]]),
            "y_pos": np.array([[0, 0, -1], [0, 1, 0], [1, 0, 0]]),
            "y_neg": np.array([[0, 0, 1], [0, 1, 0], [-1, 0, 0]]),
            "z_pos": np.array([[0, -1, 0], [1, 0, 0], [0, 0, -1]]),
            "z_neg": np.array([[0, 1, 0], [-1, 0, 0], [0, 0, -1]])
        }

    def load_data(self, command=None) -> None:
        """
        Load data into the protein library.
        If the command "manual" is not given, the input string is implemented in the data.
        Otherwise a manually folded protein is implemented in the data.
        """
        if command is None:
            for i, char in enumerate(self.sequence):
                self.data[(i, 0, 0)] = (f"{char}", i)
            return
        command.lower()
        
        if command == "manual":
            self.data[(0, 0, 0)] = ("P", 0)
            self.data[(1, 0, 0)] = ("C", 1)
            self.data[(1, 1, 0)] = ("H", 2)
            self.data[(0, 1, 0)] = ("H", 3)
            self.data[(0, 2, 0)] = ("H", 4)
            self.data[(1, 2, 0)] = ("P", 5)
            self.data[(2, 2, 0)] = ("P", 6)
            self.data[(2, 1, 0)] = ("H", 7)
            self.data[(2, 0, 0)] = ("C", 8)
            self.data[(2, -1, 0)] = ("H", 9)
            self.data[(1, -1, 0)] = ("C", 10)
            self.data[(1, -1, 1)] = ("P", 11)
            self.data = dict(sorted(self.data.items(), key=lambda item: item[1][1]))

    def give_data(self) -> dict[tuple[int, int, int], tuple[str, int]]:
        """Returns the data dictionary of a protein."""
        return self.data

    def neighbours(self, coord: tuple[int, int, int]) -> list[tuple[int, int, int]]:
        """Returns adjacent 3D coordinates of the given coordinate."""
        x_diff = [0, 1, 0, 0, -1, 0]
        y_diff = [1, 0, 0, -1, 0, 0]
        z_diff = [0, 0, 1, 0, 0, -1]
        result = [(x_diff[i] + coord[0], y_diff[i] + coord[1], z_diff[i] + coord[2]) for i in range(6)]
        return result

    def type_bond(self, coord1: tuple[int, int, int], coord2: tuple[int, int, int]) -> int:
        """
        Gives back the stability of the bond that could be formed between the two amino acids.

        Options:
        - Two "H" amino acids -> -1
        - A "H" with a "C" amino acid -> -1
        - Two "C" acids -> -5
        - Other options -> 0
        """
        if coord1 in self.data and coord2 in self.data:
            type_aminos = set()
            type_aminos.add(self.data[coord1][0])
            type_aminos.add(self.data[coord2][0])
            if type_aminos == {"H"} or type_aminos == {"C", "H"}:
                return -1
            elif type_aminos == {"C"}:
                return -5
        return 0

    def stability(self) -> int:
        """A function that calculates the stability of a protein and returns it as an integer."""
        # Filter out the "H" and "C" amino acids
        h_acids: dict[tuple[int, int, int], tuple[str, int]] = {}
        for acid in self.data.items():
            if acid[1][0] == "H" or acid[1][0] == "C":
                h_acids[acid[0]] = acid[1]

        # Loop through the "H" and "C" amino acids and look to the neighbours
        score = 0
        for acid in h_acids.items():
            friends = self.neighbours(acid[0])
            for friend in friends:
                if friend in self.data and abs(acid[1][1] - self.data[friend][1]) == 1:
                    continue
                else:
                    score += self.type_bond(acid[0], friend)
        score //= 2
        return score

    def rotate_coord(self, coord: tuple[int, int, int], pivot: tuple[int, int, int], matrix) -> tuple[int, int, int]:
        """
        Rotates the coordinate around a pivot with a rotation matrix.

        Input:
        - coord: the point that is being rotated.
        - pivot: is the point around which the coord is rotated.
        - matrix: the matrix that defines the direction of the rotation.

        Output:
        - Rotated coordinate
        """
        rel_coord = (coord[0] - pivot[0], coord[1] - pivot[1], coord[2] - pivot[2])
        v_rel_coord = np.array(rel_coord, dtype=int)
        rel_rot_coord = matrix @ v_rel_coord
        rel_new_coord = tuple(rel_rot_coord)
        new_coord = (rel_new_coord[0] + pivot[0], rel_new_coord[1] + pivot[1], rel_new_coord[2] + pivot[2])
        return new_coord

    def is_foldable(self, pivot: tuple[int, int, int], matrix) -> bool:
        """Returns if the protein can rotate in a certain direction (matrix) around a pivot point."""
        # Store the points that are following in the sequence
        points_to_rot: dict[tuple[int, int, int], tuple[str, int]] = {}
        
        for coord, amino in self.data.items():
            if amino[1] > self.data[pivot][1]:
                points_to_rot[coord] = amino

        # Check if the rotation doesn't clash with the existing folding
        for acid in points_to_rot:
            new_coord = self.rotate_coord(acid, pivot, matrix)
            if new_coord in self.data:
                return False
        return True
    
    def possible_folds(self) -> list[tuple[tuple[int, int, int], str]]:
        """Returns the possible folds in a protein."""
        possibilities: list[tuple[tuple[int, int, int], str]]= []
        for coord, amino in self.data.items():
            if amino[1] != 0 and amino[1] != len(self.sequence) - 1:
                for key, value in self.rotations.items():
                    if self.is_foldable(coord, value):
                        possibilities.append((coord, key))
        return possibilities
    
    def possible_folds_point(self, coord: tuple[int, int, int]) -> list[str]:
        possibilities: list[str] = []
        for direction in self.rotations:
            if self.is_foldable(coord, self.rotations[direction]):
                possibilities.append(direction)
        return possibilities
    
    def fold_reverse(self, pivot: tuple[int, int, int], direction: str) -> bool:
        """Folds a protein in reversee to the direction that is given"""
        if "pos" in direction:
            direction = direction.replace("pos", "neg")
        else:
            direction = direction.replace("neg", "pos")
        return self.fold(pivot=pivot, direction=direction)

    def fold(self, pivot: tuple[int, int, int], direction: str) -> bool:
        """Folds a protein at a given pivot point in a certain direction: "{axis}_{direction}"."""
        # Raise an error if the coordinate is not in the dataset
        if pivot not in self.data:
            print("coordinate not in dataset")
            raise IndexError

        # Make the command case insensitive and define the rotation matrix, otherwise return False
        direction = direction.lower()

        if direction in self.rotations:
            rot_matrix = self.rotations[direction]
        else:
            return False

        # Rotate every point
        if self.is_foldable(pivot, rot_matrix):
            # Store points that must be rotated
            points_to_rot: dict[tuple[int, int, int], tuple[str, int]] = {}
        
            for coord, amino in self.data.items():
                if amino[1] > self.data[pivot][1]:
                    points_to_rot[coord] = amino

            for acid in points_to_rot:
                new_coord = self.rotate_coord(acid, pivot, rot_matrix)
                store_value = self.data.pop(acid)
                self.data[new_coord] = store_value
            return True
        else:
            return False
        
    def fold_by_DNA(self, DNA: int, index: int) -> bool:
        if DNA == 0:
            return True

        DNA_to_direction = {
            1: "x_pos",
            2: "x_neg",
            3: "y_pos",
            4: "y_neg",
            5: "z_pos",
            6: "z_neg"
        }
        direction = DNA_to_direction[DNA]
        pivot = (0, 0, 0)        
        for key, value in self.data.items():
            if value[1] == index:
                pivot = key
        
        return self.fold(pivot=pivot, direction=direction)


    def check_direction(self, coord1: tuple[int, int, int], coord2: tuple[int, int, int]) -> int:
        """
        Function to check in what direction the protein is folded.
        """
        # If the x-coordinates and z-coordinates are the same, look at the y-coordinate
        if coord2[0] == coord1[0] and coord2[2] == coord1[2]:
            if coord2[1] - coord1[1] == 1:
                return -2
            else:
                return 2
        # Else the y-coordinates and z-coordinates are the same, look at the x-coordinate
        elif coord2[1] == coord1[1] and coord2[2] == coord1[2]:
            if coord2[0] - coord1[0] == 1:
                return -1
            else:
                return 1
        else:
            if coord2[2] - coord1[2] == 1:
                return -3
            else:
                return 3

    def output(self) -> list[list[str]]:
        """
        Store the configuration of the folded protein in the right format.
        """
        fold_commands = [["amino", "fold"]]
        previous_amino = (0, 0, 0)
        for amino in self.data:
            if amino == (0, 0, 0):
                amino_char = f"{self.data[amino][0]}"
                continue

            direction = self.check_direction(amino, previous_amino)
            fold_commands.append([amino_char, f"{direction}"])
            amino_char = f"{self.data[amino][0]}"
            previous_amino = amino

        fold_commands.append([f"{self.data[previous_amino][0]}", "0"])
        fold_commands.append(["score", f"{self.stability()}"])

        return fold_commands
