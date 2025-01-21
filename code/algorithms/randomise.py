from code.classes.protein import Protein
# from code.visualisation.visualisation import *
from code.algorithms.algorithm import Algorithm
import random

class Random_fold(Algorithm):
    def random_sequence(self) -> list[int]:
        """
        Returns a list of random integers that represent folds.
        """
        return [random.randint(0, 6) for i in range(len(self.protein.sequence) - 2)]

    def run(self):
        """
        Randomly folds a protein and returns the folded protein.
        """
        self.fold_sequence = self.random_sequence()
        return self.random_fold(self.protein)

    def random_fold(self, protein: Protein) -> Protein:
        """
        Randomly folds a protein multiple times and returns the folded protein.
        """
        sequence = protein.sequence
        copy_protein = Protein(sequence)

        # Loop over the amino acids in the protein
        for i in range(1, len(copy_protein.data) - 1):
            current_coord = (0, 0, 0)
            for coordinate, (amino, index) in copy_protein.data.items():
                if index == i:
                    current_coord = coordinate
                    break
            
            # Random fold choice
            fold_direction = self.fold_sequence[i - 1]

            if copy_protein.is_foldable(current_coord, copy_protein.rotations[fold_direction]):
                copy_protein.fold(current_coord, fold_direction)

        return copy_protein
