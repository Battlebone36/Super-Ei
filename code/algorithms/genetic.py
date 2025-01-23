from code.classes.protein import Protein
from code.algorithms.randomise import Random_fold
from code.algorithms.algorithm import Algorithm
import copy
from code.visualisation.visualisation import *
import random

class Genetic(Algorithm):
    def weighted_choice_parents(self, population: dict[int, tuple[Protein, list[int]]]) -> list[int]: 
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

    def tournament_selection(self, population: dict[int, tuple[Protein, list[int]]]):
        """
        Selects a parent using tournament selection. 
        """
        tournament_size = len(population) // 2
        tournament = random.sample(list(population.values()), k=tournament_size)
        return min(tournament, key=lambda p: p[0].stability())

    def build_offsprings(self, parent1: tuple[Protein, list[int]], parent2: tuple[Protein, list[int]]) -> tuple[list[int], list[int]]:
        """
        Create two new folded proteins (offsprings) by mixing the fold sequences of the two parents.
        """
        # Mix the two parents
        split_point = random.choice(range(len(parent1[1])))
        child1_folds = parent1[1][:split_point] + parent2[1][split_point:]
        child2_folds = parent2[1][:split_point] + parent1[1][split_point:]

        return child1_folds, child2_folds

    def fold_protein_by_sequence(self, protein: Protein, folds: list[int]) -> Protein:
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

    def mutation_on_sequence(self, folds: list[int], mutation_probability: float) -> list[int]:
        """
        Applies random mutations to a sequence of folds.
        """
        mutated_folds = folds.copy()

        # Apply random mutations
        for i in range(len(mutated_folds)):
            if random.random() < mutation_probability:
                mutated_folds[i] = random.randint(0, 6)

        return mutated_folds


    def run(self, protein: Protein) -> Protein:
        """
        Genetic algorithm that mimics natural selection to find the optimal folded protein.
        """
        print("Starting algorithm")
        population: dict[int, tuple[Protein, list[int]]] = {}
        population_size: int = 3*len(protein.sequence)
        generations: int = 1000
        mutation_probability: float = 0.02
        protein_id = 0

        # Create initial population - generation 0
        for i in range(population_size):
            # Create random folds and apply it to the protein
            random_protein = Random_fold(protein)
            folded_protein = random_protein.run()
            folds = random_protein.fold_sequence

            copy_folded_protein = copy.deepcopy(folded_protein)

            population[protein_id] = (copy_folded_protein, folds)
            protein_id += 1

        print("initial population created")

        # Track best solution
        best_protein = min(population.values(), key=lambda p: p[0].stability())
        generations_without_change = 0

        for generation in range(generations):
            print(f"Generation {generation + 1}/{generations}")

            # # Select the two best proteins from the population based on stability
            # selected_parents = weighted_choice_parents(population)
            # parent1 = population[selected_parents[0]]
            # parent2 = population[selected_parents[1]]

            # Tournament selection
            parent1 = self.tournament_selection(population)
            parent2 = self.tournament_selection(population)

            # Create the offsprings with crossover
            child1_folds, child2_folds = self.build_offsprings(parent1, parent2)
            
            # Apply mutations
            child1_folds = self.mutation_on_sequence(child1_folds, mutation_probability)
            child2_folds = self.mutation_on_sequence(child2_folds, mutation_probability)

            # Apply the folds to the offsprings
            child1 = Protein(protein.sequence)
            child2 = Protein(protein.sequence)

            folded_child1 = self.fold_protein_by_sequence(child1, child1_folds)
            folded_child2 = self.fold_protein_by_sequence(child2, child2_folds)

            population[protein_id] = (folded_child1, child1_folds)
            protein_id += 1

            population[protein_id] = (folded_child2, child2_folds)
            protein_id += 1

            # Find the worst two proteins in the population
            sorted_population = sorted(population.items(), key=lambda p: p[1][0].stability(), reverse=True)
            keys_worst_proteins = [sorted_population[0][0], sorted_population[1][0]]

            # Replace the worst proteins with the offsprings if they're better
            children = [(folded_child1, child1_folds), (folded_child2, child2_folds)]
            for i in range(len(keys_worst_proteins)):
                if children[i][0].stability() < population[keys_worst_proteins[i]][0].stability():
                    population[keys_worst_proteins[i]] = children[i]

            # Update best solution
            current_best_protein = min(population.values(), key=lambda p: p[0].stability())
            if current_best_protein[0].stability() < best_protein[0].stability():
                best_protein = current_best_protein
                generations_without_change = 0
            else:
                generations_without_change += 1

            print(f"Best protein stability at generation {generation + 1}: {best_protein[0].stability()}")

            # Early exit if there is no improvement for consecutive generations
            if generations_without_change > 100:
                print(f"No improvement for 100 generations. Early stopping at generation {generation + 1}.")
                break

        print(f"Final best protein stability: {best_protein[0].stability()}")

        self.protein = best_protein[0]
        return best_protein[0]


if __name__ == "__main__":
    test = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")
    # prot = Random_fold(test)
    # prot.run()
    gen = Genetic(test)
    gen.run()
    gen.visualise()
    # visualise_protein(best_protein)