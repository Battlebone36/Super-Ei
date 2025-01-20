from code.classes.protein import Protein
# from code.visualisation.visualisation import *
from code.algorithms.algorithm import Algorithm
import random

class Random_fold(Algorithm):
    def random_sequence(self) -> list[int]:
        """
        Returns a list of random integers that represent folds.
        """
        return [random.randint(0, 6) for i in range(len(self.protein.sequence) - 2)]

    def valid_random_sequence(self) -> list[int]:
        """
        Returns a list with a valid sequence of folds.
        """
        # Make a valid fold sequence while there is none
        folds = []
        while len(folds) != len(self.copy_protein.sequence) - 2:

            # Straighten the check protein and loop through the data
            self.copy_protein.load_data()
            for coordinate, amino in self.copy_protein.data.items():

                # Skip the first and the last amino acid
                if amino[1] == 0 or amino[1] == len(self.copy_protein.sequence) - 1:
                    continue

                # Look for possible folds and add them to the list
                possible_folds = self.copy_protein.possible_folds_point(coordinate)
                if possible_folds:
                    choice = random.choice(possible_folds)
                    folds.append(choice)

                    # Fold the check protein
                    self.copy_protein.fold(coordinate, choice)
                
                # If there is no option stop and try again
                else:
                    break
        return folds
    
    
    def fold_protein_by_sequence(self, fold_sequence: list[int]) -> Protein:
        """
        Folds a protein by a sequence of integers
        """
        # Straighten the check protein and loop through the data
        self.copy_protein 
        for i, (coordinate, amino) in enumerate(self.copy_protein.data.items()):
            
            # Skip the first and last amino acid
            if amino[1] == 0 or amino[1] == len(self.copy_protein.sequence) - 1:
                continue

            # Check if the fold is possible and fold it
            fold_direction = fold_sequence[i - 1]
            if self.copy_protein.is_foldable(coordinate, self.copy_protein.rotations[fold_direction]):
                self.copy_protein.fold(coordinate, fold_direction)
        return self.copy_protein

    def run(self):
        """
        Randomly folds a protein and returns the folded protein.
        """
        sequence = self.protein.sequence
        copy_protein = Protein(sequence)

        fold_sequence = self.valid_random_sequence()
        return self.fold_protein_by_sequence(fold_sequence=fold_sequence)


def random_fold(protein: Protein) -> Protein:
    """
    Randomly folds a protein multiple times and returns the folded protein.
    """
    sequence = protein.sequence
    copy_protein = Protein(sequence)

    # Keep folding the protein randomly
    # directions = ["x_pos", "x_neg", "y_pos", "y_neg", "z_pos", "z_neg"]
    directions = [i for i in range(7)]

    # for attempt in range(int(len(protein.data) * 1.5)):
        

    # Loop over the amino acids in the protein
    for i in range(1, len(copy_protein.data) - 1):
        current_coord = (0, 0, 0)
        for coord, (amino, index) in copy_protein.data.items():
            if index == i:
                current_coord = coord
                break
        
        # Random fold choice
        fold_direction = random.choice(directions)

        if copy_protein.is_foldable(current_coord, copy_protein.rotations[fold_direction]):
            copy_protein.fold(current_coord, fold_direction)
            
    return copy_protein

class random_fold_grow:
    def random_fold2(protein: Protein) -> Protein:
        sequence = protein.sequence
        new_protein = Protein("", "manual")
        cursor_coordinate = (0, 0, 0)
        for i, char in enumerate(sequence):
            if i == 0:
                new_protein.add_amino(cursor_coordinate, char, i)
            else:
                options = random_fold_grow.options_to_grow(new_protein, cursor_coordinate)
                if options:
                    move = random.choice(options)
                    place_coord = random_fold_grow.add_two_tuples(move, cursor_coordinate)
                    new_protein.add_amino(place_coord, char, i)
                    cursor_coordinate = place_coord
                else:
                    return protein
        return new_protein

    def add_two_tuples(tuple1: tuple[int, int, int], tuple2: tuple[int, int, int]) -> tuple[int, int, int]:
        return tuple(a + b for a, b in zip(tuple1, tuple2))

    def subtract_two_tuple(tuple1: tuple[int, int, int], tuple2: tuple[int, int, int]) -> tuple[int, int, int]:
        return tuple(a - b for a, b in zip(tuple1, tuple2))
    
    def options_to_grow(new_protein: Protein, coord: tuple[int, int, int]) -> list[tuple[int, int, int]]:
        all_options = new_protein.neighbours(coord=coord)
        filtered_options: list[tuple[int, int, int]] = []
        for option in all_options:
            if not new_protein.is_in_data(option):
                rel_coord = random_fold_grow.subtract_two_tuple(option, coord)
                filtered_options.append(rel_coord)
        return filtered_options


