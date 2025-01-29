from code.algorithms.randomise import Random
from code.classes.protein import Protein
from code.visualisation.visualisation import *
import copy

class DepthFirst(Random):
    def next_fold_sequence(self):
        """
        Advances the fold sequence by one step, similar to a clock increment.

        The fold sequence is a list of integers where each integer can range from 0 to 6.
        This method increments the sequence by one, handling overflow such that if an 
        element exceeds 6, it resets to 0 and the overflow is added to the next element 
        to the left.

        For example:
        - [0, 0, 1] will turn into [0, 0, 2]
        - [0, 0, 6] will turn into [0, 1, 0]
        - [6, 6, 6] will turn into [0, 0, 0]

        If the entire sequence is at its maximum value (all elements are 6), the next 
        sequence will reset all elements to 0.
        """

        # Add 1 to the last integer in the sequence and loop back to the beginning 
        # add the overload to the following integers
        carry = 1
        for i in reversed(range(len(self.fold_sequence))):
            self.fold_sequence[i] += carry
            carry = self.fold_sequence[i] // 7
            self.fold_sequence[i] = self.fold_sequence[i] % 7

    def index_of_invalid_fold(self) -> int:
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

    def prune_sequence(self, index: int) -> None:
        """
        Prunes the fold sequence at the specified index and adjusts the sequence accordingly.
        This method modifies the fold sequence by incrementing the value at the specified index
        and resetting subsequent values. If the increment results in a carry-over beyond the 
        maximum allowed value (6), it propagates the carry to the preceding elements. If the 
        carry reaches the beginning of the sequence, the entire fold sequence is reset to zeros.

        Args:
        ------
        index (int): The index at which to prune the fold sequence. Must be within the range 
        of the sequence length minus two.
        """
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
            self.fold_sequence = [0] * (len(self.protein.sequence) - 2)
    
    def load_possible_folds(self, verbose: bool):
        """
        Load all possible folds in a list.
        """
        # Define the first fold sequence and variables
        self.fold_sequence = [0 for i in range(len(self.protein.sequence) - 2)]
        self.next_fold_sequence()
        self.iterations = 0

        # Loop through all the possible fold sequences
        while self.fold_sequence != [0 for i in range(len(self.protein.sequence) - 2)]:
            
            # Pick the next fold seguence and prune_sequence if there is an error
            self.next_fold_sequence()
            index = self.index_of_invalid_fold()
            if index != -1:
                self.prune_sequence(index=index)
            
            # Add the fold sequence to the list
            self.storage_fold_sequences.append(tuple(self.fold_sequence))

            # Print the status if asked for
            if self.iterations % 2000 == 0 and verbose:
                print(f"Loaded {self.iterations} sequences")
            self.iterations += 1

    def evaluate_fold_sequences(self, verbose: bool = False, store_iteration_stability: bool = False) -> Protein:
        """
        Evaluates all fold sequences stored in the object and returns the protein with the best stability.

        Args:
        -------
        verbose (bool): If True, prints detailed information about the evaluation process. Default is False.
        store_iteration_stability (bool): If True, stores the stability of the protein at each iteration. Default is False.

        Returns:
        ----------
        Protein: The protein object with the best (lowest) stability found during the evaluation.
        The method iterates over all fold sequences stored in `self.storage_fold_sequences`, evaluates the stability of 
        each resulting protein, and keeps track of the protein with the best stability. If `verbose` is enabled, it prints 
        the new best stability whenever a better protein is found and periodically prints the progress percentage. If 
        `store_iteration_stability` is enabled, it stores the stability of the protein at each iteration.
        """
        """
        Evaluates all fold sequences and returns the protein with the best stability.
        """
        best_protein = self.protein
        best_stability = 0
        self.iterations = 0

        for fold_sequence in self.storage_fold_sequences:
            self.fold_sequence = fold_sequence
            candidate_protein = self.fold_by_sequence()
            candidate_stability = candidate_protein.stability()

            if candidate_stability < best_stability:
                best_protein = copy.deepcopy(candidate_protein)
                best_stability = candidate_stability
                if verbose:
                    print(f"New best stability: {best_stability}")

            if verbose and self.iterations % 2000 == 0:
                status = 0
                for i, item in enumerate(reversed(self.fold_sequence[:-3])):
                    status += item * 6 ** i
                max = 0
                for i in range(len(self.fold_sequence[:-3])):
                    max += 6 * 6 ** i
                print(f"{status * 100 // max}%")
            self.iterations += 1

            if store_iteration_stability:
                self.store_iteration_stability()

        return best_protein

    def run(self, verbose: bool = False, store_iteration_stability: bool = False) -> Protein:
        """
        Runs the Depth First Search algorithm to find the best protein fold.
        """
        self.storage_fold_sequences = []
        self.load_possible_folds(verbose=verbose)
        self.protein = self.evaluate_fold_sequences(verbose=verbose, store_iteration_stability=store_iteration_stability)
        return self.protein