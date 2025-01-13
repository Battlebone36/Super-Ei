from code.classes.protein import Protein
from code.visualisation.visualisation import *
from code.algorithms.randomise import random_fold
import csv

if __name__ == "__main__":
    
    protein_vis = Protein("HHCPPPPH", "manual")
    # visualise_protein(protein_vis)

    # plot_algorithm(random_fold)
    algorithms = [random_fold, random_fold]
    visualise_algorithm(algorithms)


    with open ('output.csv','w',newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ',')
        my_writer.writerows(protein_vis.output())
