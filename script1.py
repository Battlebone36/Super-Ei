from code.data.data_store import *


if __name__ == "__main__":
    algorithms = [Random_fold, Greedy, Climbing_fold, SimulatedAnnealing, Genetic]
    visualise_algorithm_data(algorithms)
