from code.classes.protein import Protein
from code.visualisation.visualisation import *
from code.algorithms.randomise import Random_fold
from code.algorithms.algorithm import Algorithm
import copy


class Climbing_fold(Algorithm):
    def run(self, store_step_stability: bool=False):
        """
        Start from a random state which then receives small changes
        these small changes are the neighbouring states. The best neighbouring 
        state will be the new state. Repeat these small changes until there are
        no more good changes.
        """
        # Create the random protein
        protein_random = Random_fold(self.protein)
        prot = protein_random.run()
        self.protein = copy.deepcopy(prot)
        # visualise_protein(self.protein)

        self.solve_protein(store_step_stability)
        return self.protein

    def best_move(self, max_iterations: int, store_step_stability: bool=False):
        """
        Find the best possible move at a certain configuration.
        """
        # Loop over the randomized protein
        lowest_stability = self.protein.stability()
        adjust_protein = copy.deepcopy(self.protein)
        best_protein = copy.deepcopy(self.protein)
        for amino in self.protein.data:
            # Copy the first protein and find the possible folds at a point
            possible_folds = self.protein.possible_folds_point(amino)
            if self.iterations >= max_iterations:
                break
            # Find the stability for each of the moves and save the best
            for p_folds in possible_folds:
                adjust_protein.fold(amino, p_folds)
                stability = adjust_protein.stability()
                self.iterations += 1
                if store_step_stability:
                    self.store_steps_stability()
                if self.iterations >= max_iterations:
                    break
                best_protein, lowest_stability = self.most_stable_protein(
                    stability,
                    lowest_stability,
                    adjust_protein,
                    amino,
                    p_folds,
                    best_protein)
                

        self.protein = copy.deepcopy(best_protein)

    def solve_protein(self, store_step_stability: bool=False):
        """
        Run the climbing algorithm at depth 1 2 or 3.
        """
        max_iterations = self.max_iterations
        for i in range(100):
            stability = self.protein.stability()
            self.best_move(max_iterations, store_step_stability)
            new_stability = self.protein.stability()
            if store_step_stability:

                if self.iterations >= max_iterations and stability == new_stability:
                    print(f"The top of the hill has been found at {i} climbs")
                    break
            else:
                if stability == new_stability:
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
