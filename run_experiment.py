from code.data.data_store import *
import time


if __name__ == "__main__":

    algorithms = [Random,
                  HillClimb,
                  MountainClimb,
                  SimulatedAnnealing,
                  Genetic]
    
    # -------------------------------------------
    # for algorithm in algorithms:
    #     for i in range(100):
    #         store_data(algorithm)
    # -------------------------------------------

    visualise_algorithm_data(algorithms)