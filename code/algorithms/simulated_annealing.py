from code.classes.protein import Protein
from code.algorithms.randomise import Random_fold
from code.visualisation.visualisation import *
from code.algorithms.algorithm import Algorithm
import random
import copy
import math

class SimulatedAnnealing(Algorithm):
    def run(self, store_step_stability: bool=False) -> Protein:
        """
        Starts with a randomly folded protein that will try to apply random folds.
        Each fold that improves the stability score is accepted for the next iteration.
        Sometimes worse folds are accepted, depending on the current temperature.
        """
        # Starting values 
        protein = self.protein
        initial_temp: int = 5
        cooling_rate: int = 0.99
        min_temp: int = 1
        times: int = 10

        # Track the best solution found
        random_protein = Random_fold(protein)
        current_protein = random_protein.run()
        best_protein = copy.deepcopy(current_protein)
        current_stability = current_protein.stability()
        best_stability = current_stability

        # Initialize the temperature
        current_temp = initial_temp
        
        while current_temp > min_temp and self.iterations <= 1000:
            # Try random folding 10 times per temperature
            for i in range(times):
                possible_folds = current_protein.possible_folds()

                if not possible_folds and self.iterations <= 1000:
                    break
                else:
                    self.iterations += 1
                    if store_step_stability:
                        self.store_steps_stability()

                # Select a random fold
                pivot, direction = random.choice(possible_folds)
                
                # Attempt to apply the random fold
                new_protein = copy.deepcopy(current_protein)
                
                if new_protein.is_foldable(pivot, new_protein.rotations[direction]):
                    new_protein.fold(pivot, direction)
                    new_stability = new_protein.stability()

                    # Calculate the change in stability
                    delta_e = new_stability - current_stability

                    # Calculate the acceptance probability
                    probability = math.exp(-delta_e / current_temp)

                    # Accept or deny the random fold
                    if (delta_e < 0 or random.uniform(0, 1) < probability):
                        # New configuration is accepted
                        current_protein = new_protein
                        current_stability = new_stability

                        # Best solution is updated
                        if current_stability < best_stability:
                            best_protein = copy.deepcopy(current_protein)
                            best_stability = current_stability

            # Lower the current temperature
            current_temp *= cooling_rate
        
        self.protein = best_protein
        return best_protein