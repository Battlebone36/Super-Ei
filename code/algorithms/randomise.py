from code.classes.protein import Protein
from code.visualisation.visualisation import *
import random
# random.seed(0)


def random_fold(protein: Protein) -> Protein:
    """
    Randomly folds a protein multiple times and returns the folded protein.
    """
    sequence = protein.sequence
    copy_protein = Protein(sequence)

    # Keep folding the protein randomly
    directions = ["x_pos", "x_neg", "y_pos", "y_neg", "z_pos", "z_neg"]

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


