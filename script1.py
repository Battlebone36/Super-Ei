from code.data.data_store import *


if __name__ == "__main__":
    algorithms = [Mountain_fold,
                SimulatedAnnealing,
                Gen_2
                  ]
    # for i in range(100):
    #     print(i)
    #     for algorithm in algorithms:
    #         store_data(algorithm)
    visualise_algorithm_data(algorithms)
