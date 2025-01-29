from code.classes.protein import Protein
from code.algorithms.depth import DepthFirst
from code.algorithms.genetic import Genetic
from code.algorithms.greedy import Greedy
from code.algorithms.hillclimb import HillClimb
from code.algorithms.mountainclimb import MountainClimb
from code.algorithms.simulated_annealing import SimulatedAnnealing
from code.algorithms.randomise import Random

import sys

# Check for possible algorithms
if "--Print_Possible_Algorithms" in sys.argv:
    print("Possible algorithms are: Random, DepthFirst, Greedy, HillClimb, MountainClimb, SimulatedAnnealing, Genetic")
    sys.exit(0)

# Check for proper usage
if len(sys.argv) < 2:
    print("Usage: python3 fold_by_algorithm.py <Algorithm1> [<Algorithm2> ...] [--Protein_sequence <str=HHPCHHHPCHPHHCHPH> --Print_Possible_Algorithms <bool=False>]")
    raise ValueError("No algorithm given")

# Define the protein sequence
sequence = "HHPCHHHPCHPHHCHPH"
if "--Protein_sequence" in sys.argv:
    index = sys.argv.index("--Protein_sequence")
    sequence = sys.argv[index+1]

# Turn strings into classes
algorithms = sys.argv[1:]
algorithms = [algorithm.lower() for algorithm in algorithms]
algorithm_limit = len(algorithms)

for i in range(len(algorithms)):
    if algorithms[i] == "hillclimb":
        algorithms[i] = HillClimb
    elif algorithms[i] == "mountainclimb":
        algorithms[i] = MountainClimb
    elif algorithms[i] == "simulatedannealing":
        algorithms[i] = SimulatedAnnealing
    elif algorithms[i] == "genetic":
        algorithms[i] = Genetic
    elif algorithms[i] == "greedy":
        algorithms[i] = Greedy
    elif algorithms[i] == "depthfirst":
        algorithms[i] = DepthFirst
    elif algorithms[i] == "random":
        algorithms[i] = Random
    else:
        algorithm_limit -= 1

# Filter the algorithms out of the arguments
algorithms = algorithms[:algorithm_limit]


if __name__ == "__main__":

    # Run all algorithms on the protein sequence and visualise the results
    for algorithm in algorithms:
        protein = Protein(sequence)
        temp_alg = algorithm(protein)
        temp_alg.run()
        print(f"{algorithm.__name__} is done")
        temp_alg.visualise()
        
