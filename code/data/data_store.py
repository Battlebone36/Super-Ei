
from code.visualisation.visualisation import *
import csv 
from pathlib import Path

from code.classes.protein import Protein
from code.algorithms.randomise import Random_fold
from code.algorithms.depth import DepthFirst
from code.algorithms.greedy import Greedy
from code.algorithms.hillclimb import Climbing_fold
from code.algorithms.simulated_annealing import SimulatedAnnealing
from code.algorithms. genetic import Genetic

# def gather_data(algorithm) -> list[str, int, float]:
#     """
#     2nd generation helper function

#     Test the algorithm 100 times and store the data into lists.
#     Return the data.
#     """
#     sequence = "PPCHHPPCHPPPPCHCPCHC"
#     test_protein = Protein(sequence)
#     data: list[str, int, float] = []
    

#     # Test the algorithm 100 times and store the result
#     for i in range(1000):
#         print(i)
#         start = timer()
#         result_protein: Protein = algorithm(test_protein)
#         end = timer()

#         data.append([sequence, result_protein.stability(), end - start, f"{algorithm.__name__}"])

#     return data

def store_data(algorithm, data: list[list[str, int, float]]) -> None:
    """
    Store the data into a csv file. 
    """
    fname = f"code/data/csv_data/{algorithm.__name__}.csv"
    write_mode = 'w'
     
    my_file = Path(fname)
    if my_file.is_file():
        write_mode = 'a'

    with open (fname, write_mode, newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ',')
        if write_mode == 'w':
            my_writer.writerow(["protein", "stability", "time", "algorithm"])
        my_writer.writerows(data)

def run_algorithms(algorithms, times: int=100, type_data:str ="step_stability"):
    """
    Run all algorithms an amount of times with the store the step stability command turned on.
    """
    test = Protein("PPCHHPPCHPPPPCHHHHCHHPPHHPPPPHHPPHPP")
    
    if type_data == "step_stability":
        for algorithm in algorithms:
            for i in range(times):
                print(f"{i} of {times} of {algorithm.__name__}")
                temp_class = algorithm(test)
                temp_class.run(store_step_stability=True)

    elif type_data == "occurency_stability":
        data: list[list[str, int, float]] = []

        # Test the algorithm 100 times and store the result
        for i in range(times):
            print(i)
            start = timer()
            alg: Protein = algorithm(test)
            result_protein: Protein = alg.run()
            end = timer()

            data.append([test.sequence, result_protein.stability(), end - start, f"{algorithm.__name__}"])



if __name__ == "__main__":
    algorithms = [Climbing_fold, SimulatedAnnealing, Genetic]
    run_algorithms(algorithms, 1)