from code.classes.protein import Protein
from code.visualisation.visualisation import *
from code.algorithms.randomise import random_fold
from code.algorithms.greedy import greedy_fold
from code.algorithms.greedy_search import greedy_search_sequence
from code.algorithms.hillclimb import climbing_fold, better_climbing_fold, even_better_climbing_fold
from code.algorithms.simulated_annealing import simulated_annealing
import csv 
from pathlib import Path

def gather_data(algorithm) -> list[[str, int, float]]:
    """
    2nd generation helper function

    Test the algorithm 100 times and store the data into lists.
    Return the data.
    """
    sequence = "PPCHHPPCHPPPPCHCPCHC"
    test_protein = Protein(sequence)
    data: list[str, int, float] = []
    

    # Test the algorithm 100 times and store the result
    for i in range(100):
        print(i)
        start = timer()
        result_protein: Protein = algorithm(test_protein)
        end = timer()

        data.append([sequence, result_protein.stability(), end - start])

    return data

def store_data(algorithm: str) -> None:
    """
    Store the data into a csv file. 
    """
    data = gather_data(algorithm)
    fname = f"code/data/{algorithm.__name__}.csv"
    write_mode = 'w'
        
    my_file = Path(fname)
    if my_file.is_file():
        write_mode = 'a'

    with open (fname, write_mode, newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ',')
        if write_mode == 'w':
            my_writer.writerow(["protein", "stability", "time"])
        my_writer.writerows(data)

if __name__ == "__main__":
    store_data(better_climbing_fold)
    