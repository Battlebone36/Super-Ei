from code.classes.protein import Protein
from code.visualisation.visualisation import *
import random
import copy
import math
random.seed(0)

def simulated_annealing(protein: Protein, initial_temp: int = 100, cooling_rate: int = 0.99, min_temp = 1) -> Protein:
    """"""
    # Track the best solution found
    current_protein = copy.deepcopy(protein)
    best_protein = copy.deepcopy(protein)
    current_stability = current_protein.stability()
    best_stability = current_stability
    
    # Initialize the temperature
    current_temp = initial_temp
    
    while current_temp > min_temp:
        for _ in range(10):
            print(current_temp)
            possible_folds = current_protein.possible_folds()

            if not possible_folds:
                break

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

        # Lower the temperature
        current_temp *= cooling_rate

    return best_protein


# protein = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")
# best_protein = simulated_annealing(protein)
# visualise_protein(best_protein)



# delta_e < 0 -> for the better configurations -> always accepted
# random.uniform(0, 1) < probability -> for the worse configurations -> randomly accepted or denied