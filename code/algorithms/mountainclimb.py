from code.algorithms.hillclimb import HillClimb
import copy

class MountainClimb(HillClimb):
    def best_move(self, max_iterations: int, store_iteration_stability: bool=False):
        """
        Find the best possible move for the protein configuration by evaluating 
        the stability of potential folds over two moves.
        
        This method iterates over each amino acid in the protein, evaluates all 
        possible folds at that point, and then evaluates all possible folds for 
        the next move. It keeps track of the configuration with the lowest 
        stability and updates the protein accordingly.

        Args:
        -------
        max_iterations (int): The maximum number of iterations to perform.
        store_iteration_stability (bool, optional): If True, store the stability 
        of the protein at each iteration. Defaults to False.
        """
        # Loop over the randomized protein
        lowest_stability = self.protein.stability()
        adjust_protein = copy.deepcopy(self.protein)
        best_protein = copy.deepcopy(self.protein)
        loop_protein = copy.deepcopy(self.protein)
        for amino in loop_protein.data:
            # Copy the first protein and find the possible folds at a point
            possible_folds = loop_protein.possible_folds_point(amino)
            # Find the stability for each of the moves and save the best
            for p_folds in possible_folds:
                adjust_protein.fold(amino, p_folds)
                stability = adjust_protein.stability()
                # Check another move
                # -----------------------------------------------------------------
                mid_adjust_protein = copy.deepcopy(adjust_protein)
                for mid_amino in adjust_protein.data:

                    mid_possible_folds = adjust_protein.possible_folds_point(mid_amino)

                    for mid_p_folds in mid_possible_folds:
                        mid_adjust_protein.fold(mid_amino, mid_p_folds)
                        mid_stability = mid_adjust_protein.stability()
                        best_protein, lowest_stability = self.most_stable_protein(
                            mid_stability,
                            lowest_stability,
                            mid_adjust_protein,
                            mid_amino,
                            mid_p_folds,
                            best_protein)
                        self.iterations += 1
                        if self.iterations >= max_iterations:
                            break
                        if store_iteration_stability:
                            self.protein = best_protein
                            self.store_iteration_stability()
                            
                    # -----------------------------------------------------------------
                best_protein, lowest_stability = self.most_stable_protein(
                    stability,
                    lowest_stability,
                    adjust_protein,
                    amino,
                    p_folds,
                    best_protein)
                

        self.protein = copy.deepcopy(best_protein)