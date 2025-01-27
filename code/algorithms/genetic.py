from code.classes.protein import Protein
from code.algorithms.randomise import Random
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
            random_protein = Random(self.protein)
            folded_protein = random_protein.run()
            folds = random_protein.fold_sequence

            copy_folded_protein = copy.deepcopy(folded_protein)

            init_population[protein_id] = (copy_folded_protein, folds)
            protein_id += 1
        
        return init_population, protein_id

    def tournament_selection(self, population: dict[int, tuple[Protein, list[int]]], nominator: int):
        """
        Selects a parent using tournament selection. 
        """
        # Determine the amount of proteins in the tournament
        tournament_size = len(population) * nominator // 100

        # Create the tournament pool 
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


    def run(self, store_step_stability: bool=False, population_size: int = 30, mutation_probability: float = 0.01, nominator: int = 70, verbose: bool = False) -> Protein:
        """
        Genetic algorithm that mimics natural selection to find the optimal folded protein.
        """
        if verbose:
            print("Starting algorithm")

        # Initialize variables
        population: dict[int, tuple[Protein, list[int]]] = {}
        generations: int = 1000
        iteration_limit = self.max_iterations

        # Create initial population - generation 0
        population, protein_id = self.create_initial_population(population_size)
        if verbose:
            print("initial population created")

        # Track best solution
        best_protein = min(population.values(), key=lambda p: p[0].stability())
        generations_without_change = 0
        store_max_protein = self.protein


        for generation in range(generations):
            if self.iterations >= iteration_limit:
                break
            
            if verbose:
                print(f"Generation {generation + 1}/{generations}")

            # Create new population
            new_population: dict[int, tuple[Protein, list[int]]] = {}
            
            # Fill the new population with offsprings
            while len(new_population) < population_size:
                # Select the "best" two proteins from the population based on stability
                parent1 = self.tournament_selection(population, nominator)
                parent2 = self.tournament_selection(population, nominator)
                
                self.iterations += 2
                if self.iterations >= iteration_limit:
                    break
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
                temporarily_best_protein = min(new_population.values(), key=lambda p: p[0].stability())

                # Track solution for data store
                current_best_stability = store_max_protein.stability()
                new_stability = temporarily_best_protein[0].stability()
                if current_best_stability > new_stability:
                    self.protein = copy.deepcopy(temporarily_best_protein[0])
                
            population = new_population

            # Update best solution
            current_best_protein = min(population.values(), key=lambda p: p[0].stability())
            if current_best_protein[0].stability() < best_protein[0].stability():
                best_protein = current_best_protein
                generations_without_change = 0
            else:
                generations_without_change += 1

            if verbose:
                print(f"Best protein stability at generation {generation + 1}: {best_protein[0].stability()}")

            # Early exit if there is no improvement for consecutive generations

            if generations_without_change > 100 and not store_step_stability:
                if verbose:
                    print(f"No improvement for 100 generations. Early stopping at generation {generation + 1}.")
                break
        
        if verbose:
            print(f"Final best protein stability: {best_protein[0].stability()}")

        self.protein = best_protein[0]
        return best_protein[0]


if __name__ == "__main__":
    test = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")
    gen = Genetic(test)
    gen.run()
    gen.visualise()