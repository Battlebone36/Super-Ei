from code.data.data_store import *
import time


if __name__ == "__main__":
    # algorithms = [Mountain_fold,
    #             SimulatedAnnealing,
    #             Gen_2
    #               ]
    # # for i in range(100):
    # #     print(i)
    # #     for algorithm in algorithms:
    # #         store_data(algorithm)
    # visualise_algorithm_data(algorithms)


    algorithms = [Random,
                  HillClimb,
                  Mountain_fold,
                  SimulatedAnnealing,
                  Gen_2]
    for algorithm in algorithms:
        start = time.time()
        n_runs = 0
        
        while time.time() - start < 300:
            print(f"run: {n_runs}")
            experiment_data(algorithm)
            n_runs += 1