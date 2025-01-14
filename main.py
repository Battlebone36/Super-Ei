from code.classes.protein import Protein
from code.visualisation.visualisation import *
from code.algorithms.randomise import random_fold
from code.algorithms.greedy import greedy_fold
from code.algorithms.greedy_search import greedy_search_sequence
import csv

def write_output(protein: Protein):
    with open ('output.csv','w',newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ',')
        my_writer.writerows(protein.output())


if __name__ == "__main__":
    
    # protein_vis = Protein("HHCPPPPH", "manual")
    # visualise_protein(protein_vis)

    # algorithms = [random_fold, greedy_fold, greedy_search_sequence]
    algorithms = [greedy_search_sequence]
    visualise_algorithm(algorithms)

    protein1 = Protein("HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH")



    
