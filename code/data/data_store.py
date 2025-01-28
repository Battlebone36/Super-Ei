
from code.visualisation.visualisation import *
import csv 
from pathlib import Path

from code.classes.protein import Protein
from code.algorithms.randomise import Random
from code.algorithms.depth import DepthFirst
from code.algorithms.greedy import Greedy
from code.algorithms.hillclimb import HillClimb
from code.algorithms.simulated_annealing import SimulatedAnnealing
from code.algorithms. genetic import Genetic
from code.algorithms.mountainclimb import Mountain_fold


def gather_data(algorithm) -> list[str, int, float]:
    """
    2nd generation helper function
    Test the algorithm 100 times and store the data into lists.
    Return the data.
    """
    sequence = "PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP"
    test_protein = Protein(sequence)
    data: list[str, int, float] = []

    # Test the algorithm 100 times and store the result
    for i in range(1):
        # print(f"{i} of {algorithm.__name__}")
        start = timer()
        temp_class = algorithm(test_protein)
        temp_class.run()
        result_protein = temp_class.protein
        end = timer()

        data.append([sequence, result_protein.stability(), end - start, f"{algorithm.__name__}"])

    return data

def store_data(algorithm) -> None:
    """
    Store the data into a csv file. 
    """
    data = gather_data(algorithm)
    fname = f"code/data/csv_data/{algorithm.__name__}.csv"
    write_mode = 'w'

    my_file = Path(fname)
    if my_file.is_file():
        write_mode = 'a'

    with open (fname, write_mode, newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ',')
        if write_mode == 'w':
            my_writer.writerow(["Protein", "Stability", "Time", "Algorithm"])
        my_writer.writerows(data)

def run_for_steps_and_stability(algorithms, times: int):
    """
    Run all algorithms an amount of times with the store the step stability command turned on.
    """
    test = Protein("PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP")
    for algorithm in algorithms:
        for i in range(times):
            print(f"{i + 1} of {times} of {algorithm.__name__}")
            temp_class = algorithm(test)
            temp_class.run(store_step_stability=True)

def experiment_data(algorithm) -> None:
    """
    Store the data into a csv file for the experiment. 
    """
    data = gather_data(algorithm)
    fname = f"code/data/experiment_data/{algorithm.__name__}.csv"
    write_mode = 'w'

    my_file = Path(fname)
    if my_file.is_file():
        write_mode = 'a'

    with open (fname, write_mode, newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ',')
        if write_mode == 'w':
            my_writer.writerow(["Protein", "Stability", "Time", "Algorithm"])
        my_writer.writerows(data)


if __name__ == "__main__":
    algorithms = [SimulatedAnnealing]
    run_for_steps_and_stability(algorithms, 100)