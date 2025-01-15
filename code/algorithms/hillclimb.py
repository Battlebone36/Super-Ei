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
    
    return protein

def best_move(protein: Protein) -> Protein:
    """
    Find the best possible move at a certain configuration.
    """
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
            adjust_protein.fold_reverse(amino, p_folds)
    # visualise_protein(best_protein)

    return best_protein

def better_climbing_fold(protein: Protein) -> Protein:
    """
    The same as the climbing fold but this one looks into more moves
    at the same time
    """
    # Create the random protein
    protein = random_fold(protein)
    # visualise_protein(protein)

    for i in range(20):
        stability = protein.stability()
        protein = better_best_move(protein)
        new_stability = protein.stability()
        if new_stability == stability:
            break

    return protein

def better_best_move(protein: Protein) -> Protein:
    """
    Find the best possible move at a certain configuration but look
    into more moves at once.
    """
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

            # Check another move
            # -----------------------------------------------------------------
            mid_adjust_protein = copy.deepcopy(adjust_protein)
            for mid_amino in adjust_protein.data:

                mid_possible_folds = adjust_protein.possible_folds_point(mid_amino)

                for mid_p_folds in mid_possible_folds:
                    mid_adjust_protein.fold(mid_amino, mid_p_folds)
                    mid_stability = mid_adjust_protein.stability()

                    if mid_stability < lowest_stability:
                        best_protein = copy.deepcopy(mid_adjust_protein)
                        lowest_stability = mid_stability
                    mid_adjust_protein.fold_reverse(mid_amino, mid_p_folds)


            # -----------------------------------------------------------------

            if stability < lowest_stability:
                best_protein = copy.deepcopy(adjust_protein)
                lowest_stability = stability
            adjust_protein.fold_reverse(amino, p_folds)
    # visualise_protein(best_protein)

    return best_protein
    
def even_better_climbing_fold(protein: Protein) -> Protein:
    """
    The same as the climbing fold but this one looks into more moves
    at the same time
    """
    # Create the random protein
    protein = random_fold(protein)
    # visualise_protein(protein)

    for i in range(20):
        stability = protein.stability()
        protein = even_better_best_move(protein)
        new_stability = protein.stability()
        if new_stability == stability:
            break

    return protein

def even_better_best_move(protein: Protein) -> Protein:
    """
    Find the best possible move at a certain configuration but look
    into more moves at once.
    """
    # Loop over the randomized protein
    lowest_stability = protein.stability()
    adjust_protein = copy.deepcopy(protein)
    best_protein = copy.deepcopy(protein)
    for amino in protein.data:
        print(amino)
        # Copy the first protein and find the possible folds at a point
        possible_folds = protein.possible_folds_point(amino)

        # Find the stability for each of the moves and save the best
        for p_folds in possible_folds:
            adjust_protein.fold(amino, p_folds)
            stability = adjust_protein.stability()

            # Check another move
            # -----------------------------------------------------------------
            mid_adjust_protein = copy.deepcopy(adjust_protein)
            for mid_amino in adjust_protein.data:

                mid_possible_folds = adjust_protein.possible_folds_point(mid_amino)

                for mid_p_folds in mid_possible_folds:
                    mid_adjust_protein.fold(mid_amino, mid_p_folds)
                    mid_stability = mid_adjust_protein.stability()

                    # Check another move
                    # -----------------------------------------------------------------
                    mid2_adjust_protein = copy.deepcopy(mid_adjust_protein)
                    for mid2_amino in mid_adjust_protein.data:

                        mid2_possible_folds = mid_adjust_protein.possible_folds_point(mid2_amino)

                        for mid2_p_folds in mid2_possible_folds:
                            mid2_adjust_protein.fold(mid2_amino, mid2_p_folds)
                            mid2_stability = mid2_adjust_protein.stability()

                            if mid2_stability < lowest_stability:
                                best_protein = copy.deepcopy(mid2_adjust_protein)
                                lowest_stability = mid2_stability
                            mid2_adjust_protein.fold_reverse(mid2_amino, mid2_p_folds)


                    # -----------------------------------------------------------------

                    if mid_stability < lowest_stability:
                        best_protein = copy.deepcopy(mid_adjust_protein)
                        lowest_stability = mid_stability
                    mid_adjust_protein.fold_reverse(mid_amino, mid_p_folds)


            # -----------------------------------------------------------------

            if stability < lowest_stability:
                best_protein = copy.deepcopy(adjust_protein)
                lowest_stability = stability
            adjust_protein.fold_reverse(amino, p_folds)
    # visualise_protein(best_protein)

    return best_protein


# protein1 = Protein("CPPCHPPCHPPCPPHCCPCHPPCPCHPPHPC")
# even_better_climbing_fold(protein1)