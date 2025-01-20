from code.algorithms.randomise import Random_fold
from code.classes.protein import Protein
from code.visualisation.visualisation import *

class BreadthFirst(Random_fold):

    def next_fold_sequence(self):
        """
        Adds one to the sequence of the folds like a clock.
        For instance 1, 0, 0 will turn into 2, 0, 0 and
        6, 0, 0 will turn into 0, 1, 0.
        If the fold sequence is 6, 6, 6, 6, ... the next sequence will be 0, 0, 0, 0, ....
        """
        carry = 1
        for i in range(len(self.fold_sequence)):
            self.fold_sequence[i] += carry
            carry = self.fold_sequence[i] // 7
            self.fold_sequence[i] = self.fold_sequence[i] % 7
        return self.fold_sequence

    
    def index_of_error_in_fold_sequence(self) -> int:
        """
        Gives back the index of where the impossible fold occurs or -1 if ervery fold is possible.
        """
        self.protein_check.load_data()
        for coordinate, amino in self.protein.data.items():
            # print(coordinate, amino[1] - 1, len(self.fold_sequence))
            if amino[1] == 0 or amino[1] == len(self.protein.sequence) - 1 or self.fold_sequence[amino[1] - 1] == 0:
                continue
            elif self.fold_sequence[0] == 1 or self.fold_sequence[0] == 2:
                return 0
            print(coordinate, amino[1] - 1, len(self.fold_sequence))
            if not self.protein_check.fold(coordinate, self.fold_sequence[amino[1] - 1]):
                return amino[1] - 1
        return -1

    def cut_branch(self) -> None:
        """
        Cut a branch from the sequence.
        """
        index = self.index_of_error_in_fold_sequence()
        if index == -1:
            return
        carry = 1
        for i in range(len(self.fold_sequence)):
            if i >= index:
                self.fold_sequence[i] += carry
                carry = self.fold_sequence[i] // 7
                self.fold_sequence[i] = self.fold_sequence[i] % 7 
            elif i < index:
                self.fold_sequence[i] = 0

    def run(self):
        self.fold_sequence = [0 for i in range(len(self.protein.sequence) - 2)]
        self.next_fold_sequence()
        while self.fold_sequence != [0 for i in range(len(self.protein.sequence) - 2)]:
            self.next_fold_sequence()
            self.cut_branch()
            # print(self.fold_sequence)

          
        # self.next_fold_sequence()


        # if len(self.fold_sequence) != len(self.protein.sequence) - 2:
        #     return False
        # for coordinate, amino in self.protein_check.data.items():
        #     if amino[1] == 0 or amino[1] == len(self.protein_check.sequence) - 1 or self.fold_sequence[amino[1] - 1] == 0:
        #         continue
        #     print(self.fold_sequence[amino[1] - 1])
        #     print(self.protein_check.fold(coordinate, self.fold_sequence[amino[1] - 1]))
        

        # self.next_fold_sequence()
        
        # for i in range(3):
        #     print(self.fold_sequence, i)
        #     print(self.fold_sequence_is_valid())
        #     self.next_fold_sequence()

            

