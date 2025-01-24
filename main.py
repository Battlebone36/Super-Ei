from code.classes.protein import Protein
from code.algorithms.randomise import Random_fold
from code.algorithms.depth import DepthFirst
from code.algorithms.greedy import Greedy
from code.algorithms.hillclimb import Climbing_fold
from code.algorithms.simulated_annealing import SimulatedAnnealing
from code.visualisation.visualisation import *
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
    
    visualise_algorithm_data()
    # test = Protein("HCPHCHHCH")
    # rand = DepthFirst(test)

    # rand.run(verbose=True)
    # rand.visualise()
    # visualise_protein(sim.run())
    # alg.run()
    # alg.visualise()

    # algorithms = [random_fold, greedy_fold, greedy_search_sequence]
    # visualise_algorithm_data(algorithms)




