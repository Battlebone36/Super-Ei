from code.classes.protein import Protein
from code.algorithms.randomise import Random
from code.algorithms.depth import DepthFirst
from code.algorithms.greedy import Greedy
from code.algorithms.hillclimb import HillClimb
from code.algorithms.mountainclimb import Mountain_fold
from code.algorithms.simulated_annealing import SimulatedAnnealing
from code.algorithms.genetic import Genetic
from code.algorithms.Gen import Gen_2
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
    
    # visualise_algorithm_data()
    # test = Protein("HCPHPHPHCHHHHPCCPPHPPPHPPPPCPPPHPPPHPHHHHCHPHPHPHH")
    # rand = Greedy(test)
    # rand.run(verbose=True)
    # rand.visualise()
    algorithms= [Random, HillClimb, Mountain_fold, SimulatedAnnealing, Gen_2]
    visualise_algorithm_data(algorithms, type="time_stability")





