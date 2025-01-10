from code.classes.protein import Protein
import random

def random_fold(protein: Protein, attempts: int) -> Protein:
    """Randomly folds a protein multiple times."""
    sequence = ""

    for i in range(len(protein.data)):
        for coord, (amino, index) in protein.data.items():
            if index == i:
                sequence += amino
                break

    for attempt in range(attempts):
        copy_protein = Protein(sequence)

        for i in range(1, len(copy_protein.data) - 1):
            current_coord = None
            for coord, (amino, index) in copy_protein.data.items():
                if index == i:
                    current_coord = coord
                    break

            fold_choice = random.randint(0,2)

            if fold_choice == 1:
                copy_protein.is_foldable(current_coord)
                copy_protein.fold(current_coord, "left")
            elif fold_choice == 2:
                copy_protein.is_foldable(current_coord)
                copy_protein.fold(current_coord, "right")

    return protein

protein1 = Protein("HHCPPPPH")
protein1 = random_fold(protein1)
print(protein1.give_data())