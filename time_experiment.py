from code.data.data_store import *
import time


if __name__ == "__main__":

    algorithms = [Random,
                  HillClimb,
                  MountainClimb,
                  SimulatedAnnealing,
                  Genetic]
    # ------------------------------------------------------------
    # for algorithm in algorithms:
    #     start = time.time()
    #     n_runs = 0
        
    #     while time.time() - start < 3600:
    #         print(f"run: {n_runs}")
    #         experiment_data(algorithm)
    #         n_runs += 1
    # ------------------------------------------------------------
    visualise_algorithm_data(algorithms, type="one_hour_run")
    visualise_algorithm_data(algorithms, type="time_stability")
    