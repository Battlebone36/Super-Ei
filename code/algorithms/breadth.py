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
        
    def index_of_error_in_fold_sequence(self) -> int:
        """
        Gives back the index of where the impossible fold occurs or -1 if ervery fold is possible.
        """
        self.protein_check.load_data()
        for coordinate, amino in self.protein_check.data.items():
            if amino[1] == 0 or amino[1] == len(self.protein_check.sequence) - 1 or self.fold_sequence[amino[1] - 1] == 0:
                continue
            elif self.fold_sequence[0] == 1 or self.fold_sequence[0] == 2:
                return 0
            if not self.protein_check.fold(coordinate, self.fold_sequence[amino[1] - 1]):
                return amino[1] - 1
        return -1

    # def cut_branch(self) -> None:
    #     """
    #     Cut a branch from the sequence and skip to the next branch.
    #     """
    #     index = self.index_of_error_in_fold_sequence()
    #     if index == -1:
    #         return
    #     carry = 1
    #     for i in range(len(self.fold_sequence)):
    #         if i >= index:
    #             self.fold_sequence[i] += carry
    #             carry = self.fold_sequence[i] // 7
    #             self.fold_sequence[i] = self.fold_sequence[i] % 7 
    #         elif i < index:
    #             self.fold_sequence[i] = 0
    
    # def calc_radius(self, coordinate: tuple[int, int, int]):
    #     """
    #     Calculates the radius of a coordinate to the origin.
    #     """
    #     r_2 = (coordinate[0]) ** 2 + (coordinate[1]) ** 2 + (coordinate[2]) ** 2
    #     return np.sqrt(r_2)
    

    # Maybe later
    # def hash_protein(self):
    #     """
    #     Gives a rotational invariant hash for the self.protein of the algorithm.
    #     """
    #     hash_orientation = {}
    #     for coordinate, amino in self.protein.data.items():
    #         amino_orientation = (self.calc_radius(coordinate), amino[0])
    #         print(amino_orientation)
    #         if amino_orientation not in hash_orientation:
    #             hash_orientation[amino_orientation] = 1
    #         else:
    #             hash_orientation[amino_orientation] += 1

    def run(self, shout=False):
        self.fold_sequence = [0 for i in range(len(self.protein.sequence) - 2)]
        self.next_fold_sequence()
        max_stability = 0
        max_protein = self.protein
        non_valid = set()
        while self.fold_sequence != [0 for i in range(len(self.protein.sequence) - 2)]:
            self.next_fold_sequence()
            if self.index_of_error_in_fold_sequence() != -1:
                non_valid.add(tuple(self.fold_sequence))
                continue
            temp_protein = self.fold_protein_by_sequence()
            temp_stability = temp_protein.stability()
            if temp_stability < max_stability:
                max_protein = temp_protein
                max_stability = temp_stability
                if shout == True:
                    print("New score:", max_stability)
            if (
                shout == True
                and self.fold_sequence[0] // 6 == 1
                and self.fold_sequence[1] // 6 == 1
                and self.fold_sequence[2] // 6 == 1
                and self.fold_sequence[3] // 6 == 1
            ):
                print(*self.fold_sequence)
        counter = 0

        for value in sorted(non_valid):
            print(value)
            counter += 1
            if counter > 50:
                break

        return max_protein