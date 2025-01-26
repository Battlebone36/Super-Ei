from code.algorithms.hillclimb import HillClimb
import copy

class Mountain_fold(HillClimb):
    def best_move(self, max_iterations: int, store_step_stability: bool=False):
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
                # Check another move
                # -----------------------------------------------------------------
                mid_adjust_protein = copy.deepcopy(adjust_protein)
                for mid_amino in adjust_protein.data:

                    mid_possible_folds = adjust_protein.possible_folds_point(mid_amino)

                    for mid_p_folds in mid_possible_folds:
                        mid_adjust_protein.fold(mid_amino, mid_p_folds)
                        mid_stability = mid_adjust_protein.stability()
                        self.iterations += 1
                        if store_step_stability:
                            self.store_steps_stability()
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