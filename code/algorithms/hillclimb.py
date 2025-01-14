from code.classes.protein import Protein
from code.visualisation.visualisation import *
from code.algorithms.randomise import random_fold
import copy

def climbing_fold(protein: Protein) -> Protein:
    """
    Start from a random state which then receives small changes
    these small changes are the neighbouring states. The best neighbouring 
    state will be the new state. Repeat these small changes until there are
    no more good changes.
    """
    # Create the random protein
    protein = random_fold(protein)
    # visualise_protein(protein)

    for i in range(20):
        stability = protein.stability()
        protein = best_move(protein)
        new_stability = protein.stability()
        if new_stability == stability:
            break
    

def best_move(protein: Protein) -> Protein:
    # Loop over the randomized protein
    lowest_stability = protein.stability()
    adjust_protein = copy.deepcopy(protein)
    best_protein = copy.deepcopy(protein)
    for amino in protein.data:

        # Copy the first protein and find the possible folds at a point
        possible_folds = protein.possible_folds_point(amino)

        # Find the stability for each of the moves and save the best
        for p_folds in possible_folds:
            adjust_protein.fold(amino, p_folds)
            stability = adjust_protein.stability()
            if stability < lowest_stability:
                best_protein = copy.deepcopy(adjust_protein)
                lowest_stability = stability
            adjust_protein.fold_revers(amino, p_folds)
    # visualise_protein(best_protein)

    return best_protein



# protein1 = Protein("HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH")
# climbing_fold(protein1)