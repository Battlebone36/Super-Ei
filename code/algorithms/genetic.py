from code.classes.protein import Protein
from code.algorithms.randomise import Random_fold
from code.algorithms.algorithm import Algorithm
import copy
from code.visualisation.visualisation import *
import random

class Genetic(Algorithm):
    def create_initial_population(self, population_size: int) -> dict[int, tuple[Protein, list[int]]]:
        """
        Create the initial population (generation 0) of randomly folded proteins.
        """
        init_population: dict[int, tuple[Protein, list[int]]] = {}
        protein_id = 0

        for i in range(population_size):
            # Create random folds and apply it to the protein
            random_protein = Random_fold(self.protein)
            folded_protein = random_protein.run()
            folds = random_protein.fold_sequence

            copy_folded_protein = copy.deepcopy(folded_protein)

            init_population[protein_id] = (copy_folded_protein, folds)
            protein_id += 1
        
        return init_population, protein_id

    def tournament_selection(self, population: dict[int, tuple[Protein, list[int]]]):
        """
        Selects a parent using tournament selection. 
        """
        tournament_size = len(population) // 4
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


    def run(self, store_step_stability: bool=False) -> Protein:
        """
        Genetic algorithm that mimics natural selection to find the optimal folded protein.
        """
        print("Starting algorithm")
        population: dict[int, tuple[Protein, list[int]]] = {}
        population_size: int = 3*len(self.protein.sequence)
        generations: int = 1000
        mutation_probability: float = 0.02
        # protein_id = 0

        # Create initial population - generation 0
        population, protein_id = self.create_initial_population(population_size)
        print("initial population created")

        # Track best solution
        best_protein = min(population.values(), key=lambda p: p[0].stability())
        generations_without_change = 0

        for generation in range(generations):
            print(f"Generation {generation + 1}/{generations}")
            # Create new population
            new_population = {}
            
            # Fill the new population with offsprings
            while len(new_population) < population_size:
                # Select the "best" two proteins from the population based on stability
                parent1 = self.tournament_selection(population)
                parent2 = self.tournament_selection(population)
                
                self.iterations += 2
                if store_step_stability:
                    self.store_steps_stability()

                # Create the offsprings with crossover
                child1_folds, child2_folds = self.build_offsprings(parent1, parent2)
            
                # Apply mutations
                child1_folds = self.mutation_on_sequence(child1_folds, mutation_probability)
                child2_folds = self.mutation_on_sequence(child2_folds, mutation_probability)

                # Apply the folds to the offsprings
                child1 = Protein(self.protein.sequence)
                child2 = Protein(self.protein.sequence)

                folded_child1 = self.fold_protein_by_sequence(child1, child1_folds)
                folded_child2 = self.fold_protein_by_sequence(child2, child2_folds)

                # Add offsprings to the new population
                new_population[protein_id] = (folded_child1, child1_folds)
                protein_id += 1

                new_population[protein_id] = (folded_child2, child2_folds)
                protein_id += 1

            population = new_population

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