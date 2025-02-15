from code.classes.protein import Protein
from code.visualisation.visualisation import visualise_protein
from pathlib import Path
import csv

class Algorithm:
    def __init__(self, protein: Protein):
        self.protein: Protein = Protein(protein.sequence)
        self.copy_protein: Protein = Protein(protein.sequence)
        self.fold_sequence: list[int] = []
        self.storage_fold_sequences:list[list[int]] = []
        self.stability = protein.stability()
        self.iterations = 0
        self.max_iterations = 5000

    def visualise(self):
        """
        Visualises the protein made by the algorithm.
        """
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
    
    def fold_by_sequence(self, protein: Protein | None=None) -> Protein:
        """
        Folds the protein according to a predefined fold sequence.
        This method iterates over the amino acids in the protein and folds them
        based on the specified fold sequence. If no protein is provided, it uses
        the instance's protein attribute.

        Args:
        ---------
        protein (Protein, optional): The protein to be folded. If None, the instance's protein attribute is used.
        
        Returns:
        ----------
        Protein: The folded protein.
        
        Notes:
        ---------
            - The fold sequence is assumed to be a list of directions that 
              correspond to the order of amino acids in the protein.
            - The method checks if each fold is possible using the `is_foldable` 
              method before performing the fold.
        """

        if protein is None:
            protein = self.protein

        # Loop over the amino acids in the protein
        for i in range(1, len(protein.data) - 1):
            current_coord = (0, 0, 0)
            for coordinate, (amino, index) in protein.data.items():
                if index == i:
                    current_coord = coordinate
                    break
            
            # Fold for index 
            fold_direction = self.fold_sequence[i - 1]

            if protein.is_foldable(current_coord, protein.rotations[fold_direction]):
                protein.fold(current_coord, fold_direction)

        return protein
    
    def gather_iteration_stability(self) -> list[str, int, int, str]:
        """
        Gather the sequence, stability, steps of the protein and name of the algorithm and returns it.
        """
        return [[self.protein.sequence, self.protein.stability(), self.iterations, f"{self.__class__.__name__}"]]
    
    def store_iteration_stability(self):
        """
        Stores the steps and stability data into a CSV file.
        This method gathers the iteration stability data using the 
        `gather_iteration_stability` method and writes it to a CSV file 
        located in the `code/data/iteration_stability/` directory. The 
        filename is based on the class name of the instance. If the file 
        already exists, the data is add to the file; otherwise, a new file is created 
        with headers.

        The CSV file contains the following columns:
        - Protein
        - Stability
        - Iteration
        - Algorithm

        Raises:
        --------
            IOError: If there is an issue writing to the file.
        """
        data = self.gather_iteration_stability()
        fname = f"code/data/iteration_stability/{self.__class__.__name__}.csv"
        write_mode = 'w'
        my_file = Path(fname)

        if my_file.is_file():
            write_mode = 'a'

        with open (fname, write_mode, newline = '') as csvfile:
            my_writer = csv.writer(csvfile, delimiter = ',')
            if write_mode == 'w':
                my_writer.writerow(["Protein", "Stability", "Iteration", "Algorithm"])
            my_writer.writerows(data)

    def write_output(self):
        """
        Writes the configuration of a folded protein to a CSV file in a specific format.
        """
        with open ('output.csv','w',newline = '') as csvfile:
            my_writer = csv.writer(csvfile, delimiter = ',')
            my_writer.writerows(self.protein.output())