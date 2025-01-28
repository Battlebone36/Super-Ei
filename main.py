from code.classes.protein import Protein
from code.algorithms.randomise import Random
from code.algorithms.depth import DepthFirst
from code.algorithms.greedy import Greedy
from code.algorithms.hillclimb import HillClimb
from code.algorithms.mountainclimb import Mountain_fold
from code.algorithms.simulated_annealing import SimulatedAnnealing
from code.algorithms.genetic import Genetic
from code.visualisation.visualisation import *


if __name__ == "__main__":
    algorithms= [Random, HillClimb, Mountain_fold, SimulatedAnnealing, Genetic]
    visualise_algorithm_data(algorithms, type="time_stability")
