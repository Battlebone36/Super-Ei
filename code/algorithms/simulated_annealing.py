from code.classes.protein import Protein
from code.algorithms.randomise import Random
from code.visualisation.visualisation import *
from code.algorithms.algorithm import Algorithm
import random
import copy
import math

class SimulatedAnnealing(Algorithm):
    def run(self, store_iteration_stability: bool=False) -> Protein:
        """
        Executes the simulated annealing algorithm to find an optimal protein folding configuration.
        The algorithm starts with a randomly folded protein and iteratively applies random folds.
        Occasionally, worse folds are accepted based on the current temperature to escape local minima.
        
        Args:
        ------
        store_iteration_stability (bool): If True, stores the stability score at each iteration.
        
        Returns:
        --------
        Protein: The best protein configuration found with the lowest stability score.
        
        Attributes:
        -----------
        protein (Protein): The initial protein configuration.
        max_iterations (int): The maximum number of iterations allowed.
        iterations (int): The current iteration count.
        """
        # Starting values 
        protein = self.protein
        initial_temp: float = 50.0
        cooling_rate: float = 0.99
        min_temp: float = 0.000001
        times: int = 3
        iterations_limit = self.max_iterations

        # Track the best solution found
        random_protein = Random(protein)
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
                    if store_iteration_stability:
                        self.store_iteration_stability()

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
                    try:
                        probability = math.exp(-delta_e / current_temp)
                    except OverflowError:
                        probability = 0.

                    # Accept or deny the random fold
                    if (delta_e < 0 or random.uniform(0, 1) < probability):
                        # New configuration is accepted
                        current_protein = new_protein
                        current_stability = new_stability

                        # Best solution is updated
                        if current_stability < best_stability:
                            best_protein = current_protein
                            best_stability = current_stability
                            self.protein = best_protein

            # Lower the current temperature
            current_temp *= cooling_rate
        return best_protein
    

if __name__ == "__main__":
    test = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")
    # prot = Random(test)
    # prot.run()
    gen = SimulatedAnnealing(test)
    gen.run()
    gen.visualise()
    # visualise_protein(best_protein)