from code.classes.protein import Protein
from code.algorithms.randomise import Random_fold
from code.visualisation.visualisation import *
from code.algorithms.algorithm import Algorithm
from typing import Optional
import random
import copy
import math
random.seed(2)

class SimulatedAnnealing(Algorithm):
    def run(self, store_step_stability: bool=False) -> Protein:
        """
        Starts with a randomly folded protein that will try to apply random folds.
        Each fold that improves the stability score is accepted for the next iteration.
        Sometimes worse folds are accepted, depending on the current temperature.
        """
        # Starting values 
        protein = self.protein
        initial_temp: int = 15
        cooling_rate: int = 0.99
        min_temp: int = 3
        times: int = 5
        iterations_limit = 1000

        # Track the best solution found
        random_protein = Random_fold(protein)
        current_protein = random_protein.run()
        best_protein = copy.deepcopy(current_protein)
        current_stability = current_protein.stability()
        best_stability = current_stability

        # Initialize the temperature
        current_temp = initial_temp
        
        while current_temp > min_temp and self.iterations <= iterations_limit:
            # Try random folding 10 times per temperature
            for i in range(times):
                possible_folds = current_protein.possible_folds()

                if not possible_folds or self.iterations >= iterations_limit:
                    break
                else:
                    self.iterations += 1
                    if store_step_stability:
                        self.store_steps_stability()

            # Attempt a random fold
            new_protein, new_stability = self._attempt_random_fold(current_protein)
            if new_protein and self._should_accept(current_stability, new_stability, current_temp):
                current_protein = new_protein
                current_stability = new_stability

                # Update the best protein if a new best is found
                if current_stability < self.best_stability:
                    self._update_best_solution(current_protein, current_stability)

        return current_protein, current_stability

    def _attempt_random_fold(self, protein: Protein) -> tuple[Optional[Protein], Optional[float]]:
        """
        Tries to apply a random fold to the given protein.
        Returns the new protein and its stability if successful, otherwise None.
        """
        possible_folds = protein.possible_folds()
        if not possible_folds:
            return None, None

        pivot, direction = random.choice(possible_folds)
        new_protein = copy.deepcopy(protein)
        if new_protein.is_foldable(pivot, new_protein.rotations[direction]):
            new_protein.fold(pivot, direction)
            return new_protein, new_protein.stability()

        return None, None

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
                            self.protein = best_protein

            # Lower the current temperature
            current_temp *= cooling_rate
        return best_protein