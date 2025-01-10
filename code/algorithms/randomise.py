from code.classes.protein import Protein
from code.visualisation.visualisation import *
import random


def random_fold(protein: Protein, attempts: int) -> Protein:
    """Randomly folds a protein multiple times and returns the folded protein."""
    sequence = ""

    # Stores the protein sequence
    for i in range(len(protein.data)):
        for coord, (amino, index) in protein.data.items():
            if index == i:
                sequence += amino
                break

    # Keep folding the protein randomly
    directions = ["x_pos", "x_neg", "y_pos", "y_neg", "z_pos", "z_neg"]

    for attempt in range(attempts):
        copy_protein = Protein(sequence)

        # Loop over the amino acids in the protein
        for i in range(1, len(copy_protein.data) - 1):
            current_coord = None
            for coord, (amino, index) in copy_protein.data.items():
                if index == i:
                    current_coord = coord
                    break
            
            # Random fold choice
            fold_direction = random.choice(directions)

            if copy_protein.is_foldable(current_coord, getattr(copy_protein, fold_direction)):
                copy_protein.fold(current_coord, fold_direction)
            # elif fold_choice == 1 and copy_protein.is_foldable(current_coord, copy_protein.right_turn):
            #     copy_protein.fold(current_coord, "right")

    return copy_protein


protein1 = Protein("HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH")
folded_protein = random_fold(protein1, 60)
visualise_protein(folded_protein)
print(folded_protein.give_data())
