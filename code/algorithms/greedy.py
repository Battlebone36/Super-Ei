from code.classes.protein import Protein
from code.algorithms.algorithm import Algorithm
from code.visualisation.visualisation import *
import random
import copy

class Greedy(Algorithm):
    def greedy_fold(self, store_iteration_stability: bool=False) -> tuple[tuple[int, int, int], int]:
        """
        This method attempts to fold a protein by evaluating all possible folds and 
        selecting the one that results in the lowest stability. If multiple folds 
        result in the same stability, one of them is chosen at random.

        Args:
        ------
        store_iteration_stability (bool): If True, stores the stability at each 
                                            iteration of the folding process.

        Returns:
        ---------
        tuple[tuple[int, int, int], int]: A tuple containing the best fold 
                                            (represented as a tuple of three integers) 
                                            and the corresponding stability score.
        """
        # Define the folds that are possible in this state
        folds = self.copy_protein.possible_folds()
        folds_same_score: list[tuple[tuple[int, int, int], int], int] = []
        best_fold = ((0, 0, 0), 0)

        # Loop over the possibilities
        for fold in folds:

            # Fold it, calculate stability
            self.copy_protein.fold(*fold)
            stability = self.copy_protein.stability()
            self.iterations += 1

            # Store store protein and stability if stability is lower
            if self.stability > stability:
                self.stability = stability
                best_fold = fold
                folds_same_score.clear()
            elif self.stability == stability:
                folds_same_score.append(fold)

            # Store the steps and stability if asked
            if store_iteration_stability:
                self.store_iteration_stability()
            
            # Pop the protein back to previous state
            self.copy_protein.fold_reverse(*fold)
    
        # Choose a random good fold if there was no improvement
        if best_fold == ((0, 0, 0), 0) and folds_same_score:
            best_fold = random.choice(folds_same_score)

        self.protein.fold(*best_fold)
        return best_fold

    def run(self, verbose: bool = False, store_iteration_stability: bool=False) -> Protein:
        """
        A greedy algorithm that folds the protein in a short term way to maximise stability.
        """
        # Define variables and make a random fold in the starting protein
        possible_folds = self.protein.possible_folds()
        self.protein.fold(*random.choice(possible_folds))
        self.copy_protein = copy.deepcopy(self.protein)

        # Loop through the sequence and make the best fold possible
        for i in range(20):
            best_fold = self.greedy_fold(store_iteration_stability)

            # Give the status
            if verbose and (i + 1) % 10 == 0:
                print(f"{i + 1} out of {self.max_iterations}")
            self.copy_protein.fold(*best_fold)
        return self.protein
