from code.classes.protein import Protein
from code.visualisation.visualisation import visualise_protein

class Algorithm:
    def __init__(self, protein: Protein):
        self.protein: Protein = Protein(protein.sequence)
        self.copy_protein: Protein = Protein(protein.sequence)
        self.fold_sequence: list[int] = []
        self.list_fold_sequences: list[list[int]] = []

    def visualise(self):
        """
        Visualises the protein made by the algorithm.
        """
        # print(self.protein.give_data())
        visualise_protein(self.protein)
    
    def fold_sequence_is_valid(self):
        """
        Checks if the fold sequence is valid.
        """
        self.copy_protein.load_data()
        if len(self.fold_sequence) != len(self.protein.sequence) - 2:
            return False
        for coordinate, amino in self.copy_protein.data.items():
            if amino[1] == 0 or amino[1] == len(self.copy_protein.sequence) - 1 or self.fold_sequence[amino[1] - 1] == 0:
                continue
            elif not self.copy_protein.fold(coordinate, self.fold_sequence[amino[1] - 1]):
                return False
        return True
    
    def fold_by_sequence(self, protein: Protein) -> Protein:
        """
        Randomly folds a protein.
        """
        # Loop over the amino acids in the protein
        for i in range(1, len(self.protein.data) - 1):
            current_coord = (0, 0, 0)
            for coordinate, (amino, index) in self.protein.data.items():
                if index == i:
                    current_coord = coordinate
                    break
            
            # Fold for index 
            fold_direction = self.fold_sequence[i - 1]

            if self.protein.is_foldable(current_coord, self.protein.rotations[fold_direction]):
                self.protein.fold(current_coord, fold_direction)

        return self.protein
    
    