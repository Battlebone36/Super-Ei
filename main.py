from code.classes.protein import Protein
from code.visualisation.visualisation import *
from code.algorithms.randomise import random_fold
from code.algorithms.greedy import greedy_fold
from code.algorithms.greedy_search import greedy_search_sequence
from code.algorithms.hillclimb import climbing_fold, better_climbing_fold, even_better_climbing_fold
from code.algorithms.simulated_annealing import simulated_annealing
import csv

def write_output(protein: Protein):
    """
    Writes the output of a Protein object to a CSV file in a specific format.
    """
    with open ('output.csv','w',newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ',')
        my_writer.writerows(protein.output())


if __name__ == "__main__":
    protein_vis = Protein("HCPHPCPHPCHCHPH")
    print(protein_vis.fold_by_DNA(5, 4))
    # visualise_protein(even_better_climbing_fold(protein_vis))
    visualise_protein(protein_vis)

    # algorithms = [greedy_search_sequence, climbing_fold, better_climbing_fold, even_better_climbing_fold]
    algorithms = [random_fold, greedy_search_sequence, climbing_fold]
    visualise_algorithm(algorithms, "line")
    # algorithms = [greedy_search_sequence, climbing_fold, better_climbing_fold, even_better_climbing_fold]
    # algorithms = [random_fold]
    algorithms = [random_fold, simulated_annealing]
    visualise_algorithm(algorithms)

    # protein1 = Protein("HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH")
    
