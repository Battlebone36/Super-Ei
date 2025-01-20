from code.classes.protein import Protein
from code.algorithms.randomise import Random_fold
from code.algorithms.algorithm import Algorithm
import copy
from code.visualisation.visualisation import *
import random

def weighted_choice_parents(proteins: dict[int, tuple[Protein, list[int]]]) -> list[int]: 
    """
    Selects the two "best" proteins based on their stability using a weighted choice.
    Returns the index of the proteins in the dictionary.
    """
    # total_stability = sum(protein[0].stability() for protein in proteins.values())
    
    # protein stability / total stability -> chance of getting picked during the random choice
    # The higher stability, the higher the chance of it getting picked
    total_stability = 0
    for protein in proteins.values():
        total_stability += protein[0].stability()

    weights = [protein[0].stability() / total_stability for protein in proteins.values()]

    return random.choices(list(proteins.keys()), weights=weights, k=2)

def build_offsprings(parent1: tuple[Protein, list[int]], parent2: tuple[Protein, list[int]]) -> tuple[list[int], list[int]]:
    """
    Create two new folded proteins (offsprings) by mixing the fold sequences of the two parents.
    """
    # Mix the two parents
    split_point = random.choice(range(len(parent1[1])))
    child1_folds = parent1[1][:split_point] + parent2[1][split_point:]
    child2_folds = parent2[1][:split_point] + parent1[1][split_point:]

    return child1_folds, child2_folds

def fold_protein_on_sequence(protein: Protein, folds: list[int]) -> Protein:
    """
    Folds the protein by a sequence of integers.
    """
    copy_protein = Protein(protein.sequence)
    for i in range(1, len(copy_protein.data) - 1):
        current_coord = (0, 0, 0)
        for coord, (amino, index) in copy_protein.data.items():
            if index == i:
                current_coord = coord
                break

        fold_direction = folds[i - 1]

        if copy_protein.is_foldable(current_coord, copy_protein.rotations[fold_direction]):
            copy_protein.fold(current_coord, fold_direction)
    return copy_protein 


def genetic(protein: Protein) -> None:
    """"""
    proteins_dict: dict[int, tuple[Protein, list[int]]] = {}
    unchanged_count = 0
    protein_id = 0
    generations = 200

    # Initial population - generation 0
    # Create ten randomly folded proteins
    for i in range(10):
        # Create random folds and apply it to the protein
        fold_generator = Random_fold(Protein(protein.sequence))
        folds = fold_generator.valid_random_sequence()
        folded_protein = fold_generator.fold_protein_by_sequence(folds)

        proteins_dict[protein_id] = (folded_protein, folds)
        protein_id += 1

    best_protein = min(proteins_dict.values(), key=lambda p: p[0].stability())

    # Keep running until the stability doesn't change for two consecutive generations
    while unchanged_count < 200:
        # Keep track of the best protein
        current_best_protein = min(proteins_dict.values(), key=lambda p: p[0].stability())
        current_best_stability = current_best_protein[0].stability()

        # Select the two "best" proteins based on stability
        selected_parents = weighted_choice_parents(proteins_dict)
        parent1 = proteins_dict[selected_parents[0]]
        parent2 = proteins_dict[selected_parents[1]]

        # Mix the parents to create the offsprings
        child1_folds, child2_folds = build_offsprings(parent1, parent2)
        child1 = Protein(protein.sequence)
        child2 = Protein(protein.sequence)

        # Apply the folds to the offsprings
        folded_child1 = fold_protein_on_sequence(child1, child1_folds)
        folded_child2 = fold_protein_on_sequence(child2, child2_folds)

        # Check the stability
        if current_best_protein[0].stability() < best_protein[0].stability():
            best_protein = current_best_protein
            unchanged_count = 0
        else:
            unchanged_count += 1
        

        # Add the offsprings to the population
        proteins_dict[protein_id] = (folded_child1, child1_folds)
        protein_id += 1

        proteins_dict[protein_id] = (folded_child2, child2_folds)
        protein_id += 1

        # Replace the two worse proteins in the population with the new offsprings
        # But only id the new children are better, otherwise you would make it even worse

    return best_protein[0]



    # Add stability to the already existing dictionary 
            # Split the two best proteins in two -> and add these together -> based on stability (include weighted options)
                # half 1.1 + half 2.2
                # half 2.1 + half 1.2
            # Keep it going until it can't be done anymore

    # print(proteins_dict)
    # visualise_protein(proteins_dict[0])
    # visualise_protein(proteins_dict[1])
    # visualise_protein(proteins_dict[2])
    

if __name__ == "__main__":
    test = Protein("HHPHHHPHPHHHPH")
    best_protein = genetic(test)
    visualise_protein(best_protein)