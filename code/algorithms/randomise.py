from code.classes.protein import Protein
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
        return self.fold_by_sequence(self.protein)

    
