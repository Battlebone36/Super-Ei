from code.algorithms.randomise import Random_fold
from code.classes.protein import Protein
from code.visualisation.visualisation import *
import copy

class DepthFirst(Random_fold):

    def next_fold_sequence(self):
        """
        Adds one to the fold sequence like a clock to get the next sequence.
        For instance 0, 0, 1 will turn into 0, 0, 2 and
        0, 0, 6 will turn into 0, 1, 0.
        If the fold sequence is 6, 6, 6, 6, ... the next sequence will be 0, 0, 0, 0, ....
        """
        # Add 1 to the last integer in the sequence and loop back to the beginning 
        # add the overload to the following integers
        carry = 1
        for i in reversed(range(len(self.fold_sequence))):
            self.fold_sequence[i] += carry
            carry = self.fold_sequence[i] // 7
            self.fold_sequence[i] = self.fold_sequence[i] % 7

    def index_of_error_in_fold_sequence(self) -> int:
        """
        Gives back the index of where the impossible fold occurs or
        -1 if every fold is possible.
        """
        # Straighten the check protein and loop through the data
        self.copy_protein.load_data()
        for i, (coordinate, amino) in enumerate(self.copy_protein.data.items()):
            
            # Skip the first and last amino acid
            if amino[1] == 0 or amino[1] == len(self.copy_protein.sequence) - 1 or self.fold_sequence[amino[1] - 1] == 0:
                continue

            # Check if the fold is possible otherwise return the index
            if not self.copy_protein.fold(coordinate, self.fold_sequence[i - 1]):
                return i - 1
        
        # If there is no error return -1
        return -1

    def prune(self, index: int) -> None:
        """
        Cut a branch from the sequence and skip to the next branch.
        """
        # If the index is out of the range of the list return
        if index not in range(len(self.protein.sequence) - 2):
            return
        
        # Define the one that is to be added and loop through the list in revers order
        carry = 1
        for i in reversed(range(len(self.fold_sequence))):

            # At the indexes closest to 0 the carry needs to be added
            if i <= index:
                self.fold_sequence[i] += carry
                carry = self.fold_sequence[i] // 7
                self.fold_sequence[i] = self.fold_sequence[i] % 7

            # At the indexes closest to the max everything must be reset
            else:
                self.fold_sequence[i] = 0

        # If there is an overload at the end the fold sequence is set to zero's
        if carry == 1:
            self.fold_sequence = [0 for i in range(len(self.protein.sequence) - 2)]
    
    def load_possible_fold_sequences(self, verbose: bool):
        """
        Load all possible folds in a list.
        """
        # Define the first fold sequence and variables
        self.fold_sequence = [0 for i in range(len(self.protein.sequence) - 2)]
        self.next_fold_sequence()
        iterations = 0
        max = 7 ** len(self.protein.sequence)

        # Loop through all the possible fold sequences
        while self.fold_sequence != [0 for i in range(len(self.protein.sequence) - 2)]:
            
            # Pick the next fold seguence and prune if there is an error
            self.next_fold_sequence()
            index = self.index_of_error_in_fold_sequence()
            if index != -1:
                self.prune(index=index)
            
            # Add the fold sequence to the list
            self.storage_fold_sequences.append(tuple(self.fold_sequence))

            # Print the status if asked for
            if iterations % 2000 == 0 and verbose:
                print(f"{iterations} possible folds added   max: {max}")
            iterations += 1
    
    def loop_through_fold_sequences(self, verbose: bool):
        """
        Loop through the fold sequence and store the best protein.
        """
        # Define variables
        max_protein: Protein = self.protein
        max_stability = 0
        iterations = 0

        # Make every fold and look if the stability can be improved
        for fold_sequence in self.storage_fold_sequences:
            self.fold_sequence = fold_sequence
            temp_protein = self.fold_by_sequence()
            temp_stability = temp_protein.stability()

            # If there is an improvement store it
            if temp_stability < max_stability:
                max_protein = copy.deepcopy(temp_protein)
                max_stability = temp_stability

                # Print status if asked
                if verbose == True:
                    print("New score:", max_stability)
            if iterations % 2000 == 0 and verbose:
                print(f"{iterations} / {len(self.storage_fold_sequences)} done")
            iterations +=1

        return max_protein

    def run(self, verbose: bool=False, store_step_stability: bool=False):
        """
        Runs the breadth first search.
        """
        # self.storage_fold_sequences = []
        self.load_possible_fold_sequences(verbose=verbose)
        self.protein = self.loop_through_fold_sequences(verbose=verbose)

        return self.protein