from code.classes.protein import Protein
from code.visualisation.visualisation import *
from code.algorithms.randomise import Random_fold
import copy

def climbing_fold(protein: Protein) -> Protein:
    """
    Start from a random state which then receives small changes
    these small changes are the neighbouring states. The best neighbouring 
    state will be the new state. Repeat these small changes until there are
    no more good changes.
    """
    # Create the random protein
    protein_random = Random_fold(protein)
    protein_random.run()
    # visualise_protein(protein)
    solve_protein(protein_random.protein, "depth1")
    
    return protein

def depth2_climbing_fold(protein: Protein) -> Protein:
    """
    Start from a random state and look at all the neighbouring states. 
    Then look at all the neighbouring states of those. Find the best possible
    solution untill it cannot be improved anymore.
    """
    # Create the random protein
    protein_random = Random_fold(protein)
    protein_random.run()
    # visualise_protein(protein)
    solve_protein(protein_random.protein, "depth2")

    return protein
    
def depth3_climbing_fold(protein: Protein) -> Protein:
    """
    Start from a random state and look at all the neighbouring states. 
    Then look at all the neighbouring states of those. 
    Then look at all the neighbouring states of those. Find the best possible
    solution untill it cannot be improved anymore.
    """
    # Create the random protein
    protein_random = Random_fold(protein)
    protein_random.run()
    # visualise_protein(protein)

    solve_protein(protein_random.protein, "depth3")

    return protein

def best_move(protein: Protein, depth: str) -> Protein:
    """
    Find the best possible move at a certain configuration, but looks
    into 3 moves.
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
            if depth == "depth2" or depth == "depth3":
                # Check another move
                # -----------------------------------------------------------------
                mid_adjust_protein = copy.deepcopy(adjust_protein)
                for mid_amino in adjust_protein.data:

                    mid_possible_folds = adjust_protein.possible_folds_point(mid_amino)

                    for mid_p_folds in mid_possible_folds:
                        mid_adjust_protein.fold(mid_amino, mid_p_folds)
                        mid_stability = mid_adjust_protein.stability()
                        if depth == "depth3":
                            # Check another move
                            # -----------------------------------------------------------------
                            mid2_adjust_protein = copy.deepcopy(mid_adjust_protein)
                            for mid2_amino in mid_adjust_protein.data:

                                mid2_possible_folds = mid_adjust_protein.possible_folds_point(mid2_amino)

                                for mid2_p_folds in mid2_possible_folds:
                                    mid2_adjust_protein.fold(mid2_amino, mid2_p_folds)
                                    mid2_stability = mid2_adjust_protein.stability()
                                    best_protein, lowest_stability = most_stable_protein(
                                        mid2_stability,
                                        lowest_stability,
                                        mid2_adjust_protein,
                                        mid2_amino,
                                        mid2_p_folds,
                                        best_protein)
                            # -----------------------------------------------------------------
                        best_protein, lowest_stability = most_stable_protein(
                            mid_stability,
                            lowest_stability,
                            mid_adjust_protein,
                            mid_amino,
                            mid_p_folds,
                            best_protein)
                # -----------------------------------------------------------------
            best_protein, lowest_stability = most_stable_protein(
                stability,
                lowest_stability,
                adjust_protein,
                amino,
                p_folds,
                best_protein)
    # visualise_protein(best_protein)

    return best_protein

def solve_protein(protein: Protein, depth: str) -> Protein:
    """
    Run the climbing algorithm at depth 1 2 or 3.
    """
    for i in range(20):
        stability = protein.stability()
        if depth == "depth2":
            protein = best_move(protein, depth)
        elif depth == "depth3":
            protein = best_move(protein, depth)
        else:
            protein = best_move(protein, depth)
        new_stability = protein.stability()
        if new_stability == stability:
            print(f"The top of the hill has been found at {i} climbs")
            break
        visualise_protein(protein=protein)

def most_stable_protein(stability, lowest_stability, adjust_protein, amino, p_folds, best_protein):
    """
    Check whether the given protein is the best up untill now.
    """
    if stability < lowest_stability:
        best_protein = copy.deepcopy(adjust_protein)
        lowest_stability = stability
    adjust_protein.fold_reverse(amino, p_folds)

    return (best_protein, lowest_stability)

protein1 = Protein("CPPCHPPCHPPCPPHCCPCHPPCPCHPPHPC")
# depth3_climbing_fold(protein1)