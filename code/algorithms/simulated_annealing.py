from code.classes.protein import Protein
from code.algorithms.randomise import Random_fold
from code.visualisation.visualisation import *
from code.algorithms.algorithm import Algorithm
from typing import Optional
import random
import copy
import math

# class SimulatedAnnealing(Algorithm):
#     def run(self, store_step_stability: bool=False, cooling_rate: float = 0.99, begin_temp: int= 5, end_temp: int = 1, times_fold_per_temp: int = 5) -> Protein:
#         """
#         Starts with a randomly folded protein that will try to apply random folds.
#         Each fold that improves the stability score is accepted for the next iteration.
#         Sometimes worse folds are accepted, depending on the current temperature.
#         """
#         # Starting values 
#         protein = self.protein
#         iterations_limit = 1000

#         # Track the best solution found
#         random_protein = Random_fold(protein)
#         current_protein = random_protein.run()
#         best_protein = copy.deepcopy(current_protein)
#         current_stability = current_protein.stability()
#         best_stability = current_stability

#         # Initialize the temperature
#         current_temp = begin_temp
        
#         while current_temp > end_temp and self.iterations <= iterations_limit:
#             # Try random folding 10 times per temperature
#             for i in range(times_fold_per_temp):
#                 possible_folds = current_protein.possible_folds()

#                 if not possible_folds or self.iterations >= iterations_limit:
#                     break
#                 else:
#                     self.iterations += 1
#                     if store_step_stability:
#                         self.store_steps_stability()

#                 # Select a random fold
#                 pivot, direction = random.choice(possible_folds)
                
#                 # Attempt to apply the random fold
#                 new_protein = copy.deepcopy(current_protein)
                
#                 if new_protein.is_foldable(pivot, new_protein.rotations[direction]):
#                     new_protein.fold(pivot, direction)
#                     new_stability = new_protein.stability()

#                     # Calculate the change in stability
#                     delta_e = new_stability - current_stability

#                     # Calculate the acceptance probability
#                     probability = math.exp(-delta_e / current_temp)

#                     # Accept or deny the random fold
#                     if (delta_e < 0 or random.uniform(0, 1) < probability):
#                         # New configuration is accepted
#                         current_protein = new_protein
#                         current_stability = new_stability

#                         # Best solution is updated
#                         if current_stability < best_stability:
#                             best_protein = copy.deepcopy(current_protein)
#                             best_stability = current_stability
#                             self.protein = best_protein

#             # Lower the current temperature
#             current_temp *= cooling_rate
#         return best_protein

class SimulatedAnnealing(Algorithm):
    def __init__(self, protein: Protein):
        super().__init__(protein)
        self.best_protein = None
        self.best_stability = float('inf')

    def run(self, 
            store_step_stability: bool = False, 
            cooling_rate: float = 0.99, 
            start_temp: float = 5.0, 
            end_temp: float = 1.0, 
            folds_per_temp: int = 5, 
            max_iterations: int = 1000) -> Protein:
        """
        Executes the simulated annealing algorithm to optimize protein folding.
        """
        # Initialize variables
        self.best_protein = Random_fold(self.protein).run()
        self.best_stability = self.best_protein.stability()
        current_protein = copy.deepcopy(self.best_protein)
        current_stability = self.best_stability
        current_temp = start_temp

        while current_temp > end_temp and self.iterations < max_iterations:
            current_protein, current_stability = self._perform_temperature_cycle(
                current_protein, 
                current_stability, 
                current_temp, 
                folds_per_temp, 
                store_step_stability, 
                max_iterations
            )
            current_temp *= cooling_rate

        return self.best_protein

    def _perform_temperature_cycle(self, 
                                   current_protein: Protein, 
                                   current_stability: float, 
                                   current_temp: float, 
                                   folds_per_temp: int, 
                                   store_step_stability: bool, 
                                   max_iterations: int) -> tuple[Protein, float]:
        """
        Executes one temperature cycle, attempting multiple random folds.
        """
        for _ in range(folds_per_temp):
            if self.iterations >= max_iterations:
                break

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

    def _should_accept(self, current_stability: float, new_stability: float, temperature: float) -> bool:
        """
        Determines whether a new fold configuration should be accepted based on stability and temperature.
        """
        delta_e = new_stability - current_stability
        if delta_e < 0:
            return True
        return random.uniform(0, 1) < math.exp(-delta_e / temperature)

    def _update_best_solution(self, protein: Protein, stability: float) -> None:
        """
        Updates the best protein and stability found so far.
        """
        self.best_protein = copy.deepcopy(protein)
        self.best_stability = stability
        self.protein = self.best_protein
