from code.classes.protein import Protein
from code.algorithms.randomise import Random_fold
from code.visualisation.visualisation import *
import random
import copy

class Greedy(Random_fold):

    def greedy_fold(self) -> tuple[tuple[int, int, int], int]:
        # Define the folds that are possible in this state
        folds = self.copy_protein.possible_folds()
        best_fold = folds[0]

        # Loop over the possibilities
        for fold in folds:

            # Fold it, calculate stability
            self.copy_protein.fold(*fold)
            stability = self.copy_protein.stability()

            # Store store protein and stability if stability is lower
            if self.stability > stability:
                self.stability = stability
                best_fold = fold
            
            # Pop the protein back to previous state
            self.copy_protein.fold_reverse(*fold)
        self.protein.fold(*best_fold)
        return best_fold

    def run(self) -> Protein:
        """
        A greedy algorithm that folds the protein in a short term way to maximise stability.
        """
        # Make a random protein to explore the different stabilities in the loop
        random_protein = Random_fold(self.protein)
        self.protein = random_protein.run()
        self.copy_protein = copy.deepcopy(self.protein)
        old_stability = 0

        # Loop through the sequence and store the best folded protein
        for i in range(len(self.protein.sequence)):
            best_fold = self.greedy_fold()

            # Return the protein if there is no better solution
            if old_stability == self.stability:
                return self.protein

            # Define variables to start loop again in the new protein
            self.copy_protein.fold(*best_fold)
            old_stability = self.stability

        return self.protein
