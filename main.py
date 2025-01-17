from code.classes.protein import Protein
from code.visualisation.visualisation import *
from code.algorithms.randomise import random_fold, random_fold_grow
from code.algorithms.greedy import greedy_fold
from code.algorithms.greedy_search import greedy_search_sequence
from code.algorithms.hillclimb import climbing_fold, depth2_climbing_fold, climbing_fold
from code.algorithms.simulated_annealing import simulated_annealing
import csv

def write_output(protein: Protein):
    """
    Writes the configuration of a folded protein to a CSV file in a specific format.
    """
    with open ('output.csv','w',newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ',')
        my_writer.writerows(protein.output())


if __name__ == "__main__":
    # protein_vis = Protein("HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH")
    
    # visualise_protein(simulated_annealing(protein_vis))

    # greedy_search_sequence(protein_vis)
    # protein_vis.add_amino((0, 0, 0), "H", 0)
    # new_prot = random_fold_grow.algorithm(protein_vis)
    # print(protein_vis())
    # visualise_protein(random_fold_grow.algorithm(protein_vis))
    # print(protein_vis.fold_by_DNA(5, 4))
    # # visualise_protein(even_better_climbing_fold(protein_vis))
    # visualise_protein(protein_vis)

    # algorithms = [greedy_search_sequence, climbing_fold, better_climbing_fold, even_better_climbing_fold]
    # algorithms = [random_fold, greedy_search_sequence, climbing_fold]
    # visualise_algorithm(algorithms, "together")
    # algorithms = [greedy_search_sequence, climbing_fold, better_climbing_fold, even_better_climbing_fold]
    # algorithms = [random_fold, random_fold_grow.random_fold2]
    algorithms = [random_fold, greedy_fold, greedy_search_sequence]
    # visualise_algorithm(algorithms)
    visualise_algorithm_data(algorithms, "occurency-stability", "dodge")

    # protein1 = Protein("HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH")
    
