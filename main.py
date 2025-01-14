from code.classes.protein import Protein
from code.visualisation.visualisation import *
from code.algorithms.randomise import random_fold
from code.algorithms.greedy import greedy_fold
from code.algorithms.greedy_search import greedy_search_sequence
from code.algorithms.hillclimb import climbing_fold
import csv

if __name__ == "__main__":
    
    # protein_vis = Protein("HHCPPPPH", "manual")
    # visualise_protein(protein_vis)

    plot_algorithm(random_fold)
    algorithms = [random_fold, greedy_search_sequence, climbing_fold]
    visualise_algorithm(algorithms)

    # protein1 = Protein("HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH")
    # protein1 = greedy_search_sequence(protein1)
    # visualise_protein(protein1)

    # protein1 = random_fold(protein1)

    # visualise_protein(protein1)
    # protein1 = greedy_fold(protein1)
    # visualise_protein(protein1)



    with open ('output.csv','w',newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ',')
        my_writer.writerows(protein1.output())
