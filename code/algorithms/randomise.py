from code.classes.protein import Protein
from code.visualisation.visualisation import *
import random
# random.seed(0)


def random_fold(protein: Protein) -> Protein:
    """
    Randomly folds a protein multiple times and returns the folded protein.
    """
    sequence = protein.sequence

    # Keep folding the protein randomly
    directions = ["x_pos", "x_neg", "y_pos", "y_neg", "z_pos", "z_neg"]

    for attempt in range(int(len(protein.data) * 1.5)):
        copy_protein = Protein(sequence)

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