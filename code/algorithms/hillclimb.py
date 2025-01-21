from code.classes.protein import Protein
from code.visualisation.visualisation import *
from code.algorithms.randomise import Random_fold
from code.algorithms.algorithm import Algorithm
import copy


class Climbing_fold(Algorithm):
    def run(self, depth: int|None=None):
        """
        Start from a random state which then receives small changes
        these small changes are the neighbouring states. The best neighbouring 
        state will be the new state. Repeat these small changes until there are
        no more good changes.
        """
        # Create the random protein
        protein_random = Random_fold(self.protein)
<<<<<<< HEAD
        print("3.2")
        protein_random.run()
        print("3.3")
        self.protein = copy.deepcopy(protein_random.protein)
        visualise_protein(self.protein)
=======
        prot = protein_random.run()
        self.protein = copy.deepcopy(prot)
        # visualise_protein(self.protein)

>>>>>>> 19d440f8ff00bb40cef39628586b5b5ecd67feee
        self.solve_protein(depth)
        return self.protein

    def best_move(self, depth: str):
        """
        Find the best possible move at a certain configuration, but looks
        into 3 moves.
        """
        # Loop over the randomized protein
        lowest_stability = self.protein.stability()
        adjust_protein = copy.deepcopy(self.protein)
        best_protein = copy.deepcopy(self.protein)
        for amino in self.protein.data:
            # Copy the first protein and find the possible folds at a point
            possible_folds = self.protein.possible_folds_point(amino)

            # Find the stability for each of the moves and save the best
            for p_folds in possible_folds:
                adjust_protein.fold(amino, p_folds)
                stability = adjust_protein.stability()
                if depth == 2 or depth == 3:
                    # Check another move
                    # -----------------------------------------------------------------
                    mid_adjust_protein = copy.deepcopy(adjust_protein)
                    for mid_amino in adjust_protein.data:

                        mid_possible_folds = adjust_protein.possible_folds_point(mid_amino)

                        for mid_p_folds in mid_possible_folds:
                            mid_adjust_protein.fold(mid_amino, mid_p_folds)
                            mid_stability = mid_adjust_protein.stability()
                            if depth == 3:
                                # Check another move
                                # -----------------------------------------------------------------
                                mid2_adjust_protein = copy.deepcopy(mid_adjust_protein)
                                for mid2_amino in mid_adjust_protein.data:

                                    mid2_possible_folds = mid_adjust_protein.possible_folds_point(mid2_amino)

                                    for mid2_p_folds in mid2_possible_folds:
                                        mid2_adjust_protein.fold(mid2_amino, mid2_p_folds)
                                        mid2_stability = mid2_adjust_protein.stability()
                                        best_protein, lowest_stability = self.most_stable_protein(
                                            mid2_stability,
                                            lowest_stability,
                                            mid2_adjust_protein,
                                            mid2_amino,
                                            mid2_p_folds,
                                            best_protein)
                                # -----------------------------------------------------------------
                            best_protein, lowest_stability = self.most_stable_protein(
                                mid_stability,
                                lowest_stability,
                                mid_adjust_protein,
                                mid_amino,
                                mid_p_folds,
                                best_protein)
                    # -----------------------------------------------------------------
                best_protein, lowest_stability = self.most_stable_protein(
                    stability,
                    lowest_stability,
                    adjust_protein,
                    amino,
                    p_folds,
                    best_protein)

        self.protein = copy.deepcopy(best_protein)

    def solve_protein(self, depth: int):
        """
        Run the climbing algorithm at depth 1 2 or 3.
        """
        for i in range(20):
            stability = self.protein.stability()
            self.best_move(depth)
            visualise_protein(self.protein)
            new_stability = self.protein.stability()
            if new_stability == stability:
                print(f"The top of the hill has been found at {i} climbs")
                break

    def most_stable_protein(self, stability, lowest_stability, adjust_protein, amino, p_folds, best_protein):
        """
        Check whether the given protein is the best up untill now.
        """
        if stability < lowest_stability:
            best_protein = copy.deepcopy(adjust_protein)
            lowest_stability = stability
        adjust_protein.fold_reverse(amino, p_folds)

        return (best_protein, lowest_stability)

# if __name__ == "__main__":
#     protein1 = Protein("HPHPPHHPHPPHPHHPPHPH")
#     climb = Climbing_fold(protein1)
#     climb.run(3)
