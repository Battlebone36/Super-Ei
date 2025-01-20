from code.classes.protein import Protein
from code.visualisation.visualisation import *

class Algorithm:
    def __init__(self, protein: Protein):
        self.protein: Protein = Protein(protein.sequence)
        self.protein_check: Protein = Protein(protein.sequence)
        self.fold_sequence = []

    def run(self):
        self.visualise()

    def visualise(self):
        visualise_protein(self.protein)
    
    def fold_sequence_is_valid(self):
        """
        Checks if the foldsequence is doable.
        """
        self.protein_check.load_data()
        if len(self.fold_sequence) != len(self.protein.sequence) - 2:
            return False
        for coordinate, amino in self.protein_check.data.items():
            if amino[1] == 0 or amino[1] == len(self.protein_check.sequence) - 1 or self.fold_sequence[amino[1] - 1] == 0:
                continue
            elif not self.protein_check.fold(coordinate, self.fold_sequence[amino[1] - 1]):
                return False
        return True