from code.data.data_store import *


if __name__ == "__main__":
    algorithms = [Random_fold, Greedy, Climbing_fold, Mountain_fold]

    for algorithm in algorithms:
        store_data(algorithm)
    # visualise_algorithm_data(algorithms)
