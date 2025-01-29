from code.data.data_store import *
import time
import sys

if len(sys.argv) < 2:
    print("Usage: python3 time_experiment.py <Algorithm1> [<Algorithm2> ...]")
    print("       [--Print_Possible_Algorithms <bool=False> --Type_plot <str=one_hour_run>]")
    print("Example: python3 time_experiment.py Random HillClimb MountainClimb SimulatedAnnealing Genetic --Run_experiment True --Time 300")
    sys.exit(1)

if "--Print_Possible_Algorithms" in sys.argv:
    print("Possible algorithms are: Random, HillClimb, MountainClimb, SimulatedAnnealing, Genetic")
    sys.exit(0)

run_experiment = False
time_end = 100
possible_algorithms = ["random", "greedy", "hillclimb", "mountainclimb", "simulatedannealing", "genetic"]
limit_to_algorithms = 1
type_plot = "one_hour_run"

# Check if there is a Type_plot parameter with a string
if "--Type_plot" in sys.argv:
    index = sys.argv.index("--Type_plot")
    if type(sys.argv[index+1]) != type("string"):
        raise ValueError("Type_plot parameter must be a string")
    type_plot = sys.argv[index+1]
    if type_plot != "one_hour_run" and type_plot != "time":
        raise ValueError("Type_plot parameter must be 'one_hour_run' or 'time'")

# Check if all algorithms are valid
for argument in sys.argv:
    if argument.lower() in possible_algorithms:
        limit_to_algorithms += 1

# Check if there are algorithms
if limit_to_algorithms == 1:
    raise ValueError("There are no algorithms to run")

# Check if all algorithms are placed after the script name
algorithms = sys.argv[1:limit_to_algorithms]
for algorithm in algorithms:
    if sys.argv.index(algorithm) < 1 and sys.argv.index(algorithm) >= limit_to_algorithms:
        raise ValueError("All algorithms must be placed after the script name")
algorithms = [algorithm.lower() for algorithm in algorithms]

# Turn strings into classes
for i in range(len(algorithms)):
    if algorithms[i] == "random":
        algorithms[i] = Random
    elif algorithms[i] == "greedy":
        algorithms[i] = Greedy
    if algorithms[i] == "hillclimb":
        algorithms[i] = HillClimb
    elif algorithms[i] == "mountainclimb":
        algorithms[i] = MountainClimb
    elif algorithms[i] == "simulatedannealing":
        algorithms[i] = SimulatedAnnealing
    elif algorithms[i] == "genetic":
        algorithms[i] = Genetic


if __name__ == "__main__":
    if type_plot == "one_hour_run" or type_plot == "time":
        visualise_algorithm_data(algorithms, type=type_plot)
    