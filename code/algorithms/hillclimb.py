from code.classes.protein import Protein
from code.visualisation.visualisation import *
from code.algorithms.randomise import Random
from code.algorithms.algorithm import Algorithm
import copy

class HillClimb(Algorithm):
    def run(self, store_iteration_stability: bool=False):
        """
        Start from a random state which then receives small changes
        these small changes are the neighbouring states. The best neighbouring 
        state will be the new state. Repeat these small changes until there are
        no more good changes.
        """
        # Create the random protein
        protein_random = Random(self.protein)
        prot = protein_random.run()
        self.protein = copy.deepcopy(prot)

        self.solve_protein(store_iteration_stability)
        return self.protein

    def best_move(self, max_iterations: int, store_iteration_stability: bool=False):
        """
        Find the best possible move for the protein configuration using a hill climbing algorithm.

        Parameters:
        -----------
        max_iterations (int): The maximum number of iterations to perform.
        store_iteration_stability (bool): If True, store the stability at each iteration. Default is False.

        The method iterates over the amino acids in the protein and evaluates possible folds at each point.
        It calculates the stability for each possible fold and keeps track of the configuration with the lowest stability.
        The process stops if the maximum number of iterations is reached.
        """
        # Loop over the randomized protein
        lowest_stability = self.protein.stability()
        adjust_protein = copy.deepcopy(self.protein)
        best_protein = copy.deepcopy(self.protein)
        loop_protein = copy.deepcopy(self.protein)
        for amino in loop_protein.data:
            # Copy the first protein and find the possible folds at a point
            possible_folds = loop_protein.possible_folds_point(amino)
            if self.iterations >= max_iterations:
                break
            # Find the stability for each of the moves and save the best
            for p_folds in possible_folds:
                adjust_protein.fold(amino, p_folds)
                stability = adjust_protein.stability()
                self.iterations += 1
                if store_iteration_stability:
                    self.store_iteration_stability()
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

    def solve_protein(self, store_iteration_stability: bool=False):
        """
        Executes the hill climbing algorithm to optimize the protein structure.

        Parameters:
        store_iteration_stability (bool): If True, the algorithm will stop if the maximum 
                          number of iterations is reached and the stability 
                          value remains unchanged.

        The algorithm runs for a minimum of 100 iterations or until the stability value 
        remains unchanged for two consecutive iterations. If `store_iteration_stability` 
        is set to True, the algorithm will also stop if the maximum number of iterations 
        is reached and the stability value remains unchanged.
        
        During each iteration, the algorithm:
        1. Calculates the current stability of the protein.
        2. Executes the best move to potentially improve stability.
        3. Recalculates the stability of the protein.
        
        If the stability value remains unchanged for two consecutive iterations, 
        the algorithm concludes that the top of the hill has been found and stops.
        Prints a message indicating the number of climbs taken to find the top of the hill.
        """
        max_iterations = self.max_iterations
        # Run the program for at least 100 times or if the same stability value has been
        # found twice.
        for i in range(100):
            stability = self.protein.stability()
            self.best_move(max_iterations, store_iteration_stability)
            new_stability = self.protein.stability()
            # If store iteration stability is true it will stop if the max iteration
            # has been reached.
            if store_iteration_stability:

                if self.iterations >= max_iterations and stability == new_stability:
                    print(f"The top of the hill has been found at {i} climbs")
                    break
            else:
                if stability == new_stability:
                    print(f"The top of the hill has been found at {i} climbs")
                    break

    def most_stable_protein(self, stability, lowest_stability, adjust_protein, amino, p_folds, best_protein):
        """
        Check whether the given protein is the best up until now.
        """
        if stability < lowest_stability:
            best_protein = copy.deepcopy(adjust_protein)
            lowest_stability = stability
        adjust_protein.fold_reverse(amino, p_folds)

        return (best_protein, lowest_stability)
