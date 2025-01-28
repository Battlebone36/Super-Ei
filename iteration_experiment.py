from code.data.data_store import *


if __name__ == "__main__":

    algorithms = [HillClimb,
                  MountainClimb,
                  SimulatedAnnealing,
                  Genetic]
    
    # -------------------------------------------
    # for algorithm in algorithms:
    #     run_for_iterations_and_stability(algorithm, 5000)
    # -------------------------------------------

    visualise_algorithm_data(algorithms,type="iteration_stability")