from code.classes.protein import Protein
from code.visualisation.visualisation import *
from code.algorithms.randomise import Random_fold
from code.algorithms.breadth import BreadthFirst
from code.algorithms.greedy import greedy
from code.algorithms.hillclimb import Climbing_fold
# from code.algorithms.simulated_annealing import simulated_annealing
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
    test = Protein("HHPPCHHPHHHPH")
    climb = Climbing_fold(test)
    plot = climb.run()

    visualise_protein(plot)


