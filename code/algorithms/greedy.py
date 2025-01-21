from code.classes.protein import Protein
from code.algorithms.algorithm import Algorithm
from code.visualisation.visualisation import *
import random
import copy

# class greedy(Algorithm):
#     # run
#     # make greedy change
#     def run(self):



#     def greedy_choice(self) -> bool:

    

def greedy(protein: Protein) -> Protein:
    """
    A greedy algorithm that folds the protein in a short term way to maximise stability.
    """
    # Define the variabels used
    old_score = 0
    max_score = 0


    # Make a protein to explore the different stabilities in the loop
    old_protein = copy.deepcopy(protein)
    old_protein.fold(*random.choice(protein.possible_folds()))
    new_protein = old_protein
    possible_choices_counter = 0

    # Loop through the sequence and store the best folded protein
    for i in range(len(protein.sequence)):

        # Define the folds that are possible in this state and
        # a storage for folds that have the same stability
        folds = old_protein.possible_folds()
        possible_choices = []

        # Loop over the possibilities
        for fold in folds:

            # Fold it, calculate stability
            old_protein.fold(*fold)
            stab = old_protein.stability()

            # Store store protein and stability if stability is lower
            if max_score > stab:
                max_score = stab
                new_protein = copy.deepcopy(old_protein)
            
            # store a fold if it has the same result
            elif max_score == stab:
                possible_choices.append((stab, fold))
            
            # Pop the protein back to previous state
            old_protein.fold_reverse(*fold)

        # Fold if there are options with same result but not with better results otherwise return the protein
        if old_score == max_score and possible_choices:
            old_protein.fold(*possible_choices[0][1])
            new_protein = copy.deepcopy(old_protein)
            possible_choices_counter += 1
        if (old_score == max_score and not possible_choices) or possible_choices_counter == 3:
            return new_protein

        # Define variables to start loop again in the new protein
        old_protein = copy.deepcopy(new_protein)
        old_score = max_score

    return new_protein
