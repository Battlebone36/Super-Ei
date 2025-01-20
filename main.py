from code.classes.protein import Protein
from code.visualisation.visualisation import *
from code.algorithms.randomise import random_fold, random_fold_grow, Random_fold
from code.algorithms.breadth import BreadthFirst
from code.algorithms.greedy import greedy_fold
from code.algorithms.greedy_search import greedy_search_sequence
from code.algorithms.hillclimb import climbing_fold, depth2_climbing_fold
from code.algorithms.simulated_annealing import simulated_annealing
import random
random.seed(0)
import csv

def write_output(protein: Protein):
    """
    Writes the configuration of a folded protein to a CSV file in a specific format.
    """
    with open ('output.csv','w',newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ',')
        my_writer.writerows(protein.output())


if __name__ == "__main__":
    # algorithms = [random_fold, greedy_fold, greedy_search_sequence]
    # visualise_algorithm_data(algorithms, "occurency-stability") 
    # test = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")
    test = Protein("HCHHCH")
    prot = BreadthFirst(test)
    plot = prot.run(shout=True)

    visualise_protein(plot)

    # prot = climbing_fold(test)
    # visualise_protein(prot)
    # iets = BreadthFirst(test)
    # temp_prot = iets.run()
    # visualise_protein(temp_prot)


