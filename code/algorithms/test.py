from code.classes.protein import Protein
from code.algorithms.randomise import Random_fold
from code.algorithms.algorithm import Algorithm
import copy
from code.visualisation.visualisation import *
import random

def weighted_choice_parents(population: dict[int, tuple[Protein, list[int]]]) -> list[int]: 
    """
    Selects the two "best" proteins based on their stability using a weighted choice.
    Returns the index of the proteins in the dictionary.
    """
    # total_stability = sum(protein[0].stability() for protein in proteins.values())
    
    # protein stability / total stability -> chance of getting picked during the random choice
    # The higher stability, the higher the chance of it getting picked
    total_stability = 0
    for protein in population.values():
        total_stability += protein[0].stability()

    weights = [protein[0].stability() / total_stability for protein in population.values()]

    return random.choices(list(population.keys()), weights=weights, k=2)

def tournament_selection(population: dict[int, tuple[Protein, list[int]]]):
    """
    Selects a parent using tournament selection. 
    """
    tournament_size = len(population) // 2
    tournament = random.sample(list(population.values()), k=tournament_size)
    return min(tournament, key=lambda p: p[0].stability())

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

def mutation_on_sequence(folds: list[int], mutation_probability: float) -> list[int]:
    """
    Applies random mutations to a sequence of folds.
    """
    mutated_folds = folds.copy()

    # Apply random mutations
    for i in range(len(mutated_folds)):
        if random.random() < mutation_probability:
            mutated_folds[i] = random.randint(0, 6)

    return mutated_folds 


def genetic(protein: Protein) -> Protein:
    """"""
    population: dict[int, tuple[Protein, list[int]]] = {}
    unchanged_count = 0
    protein_id = 0
    mutation_probability = 0.01

    # Create initial population - generation 0
    for i in range(3*len(protein.sequence)):
        # Create random folds and apply it to the protein
        random_protein = Random_fold(protein)
        folded_protein = random_protein.run()
        folds = random_protein.fold_sequence

        copy_folded_protein = copy.deepcopy(folded_protein)

        population[protein_id] = (copy_folded_protein, folds)
        protein_id += 1

    # Track the best solution
    best_protein = min(population.values(), key=lambda p: p[0].stability())

    # (Keep running until the stability doesn't change for two consecutive generations)
    while unchanged_count < 200:
        # Keep track of the best protein
        current_best_protein = min(population.values(), key=lambda p: p[0].stability())

        # # Select the two "best" proteins based on stability
        # selected_parents = weighted_choice_parents(population)
        # parent1 = population[selected_parents[0]]
        # parent2 = population[selected_parents[1]]

        # Tournament selection
        parent1 = tournament_selection(population)
        parent2 = tournament_selection(population)

        # Mix the parents to create the offsprings
        child1_folds, child2_folds = build_offsprings(parent1, parent2)

        # Apply mutations
        child1_folds = mutation_on_sequence(child1_folds, mutation_probability)
        child2_folds = mutation_on_sequence(child2_folds, mutation_probability)

        # Apply the folds to the offsprings
        child1 = Protein(protein.sequence)
        child2 = Protein(protein.sequence)
        folded_child1 = fold_protein_on_sequence(child1, child1_folds)
        folded_child2 = fold_protein_on_sequence(child2, child2_folds)

        # Check the stability
        if current_best_protein[0].stability() < best_protein[0].stability():
            best_protein = current_best_protein
            unchanged_count = 0
        else:
            unchanged_count += 1
        
        # Add the offsprings to the population
        population[protein_id] = (folded_child1, child1_folds)
        protein_id += 1

        population[protein_id] = (folded_child2, child2_folds)
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
    test = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")
    best_protein = genetic(test)
    visualise_protein(best_protein)