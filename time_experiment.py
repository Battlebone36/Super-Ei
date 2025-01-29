from code.data.data_store import *
import time
import sys

if len(sys.argv) < 2:
    print("Usage: python3 time_experiment.py <Algorithm1> [<Algorithm2> ...]")
    print("       [--Run_experiment <bool=False> --Time <int=100> --Print_Possible_Algorithms <bool=False> --Type_plot <str=one_hour_run>]")
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

# Check if there is a Run parameter with a boolean
if "--Run_experiment" in sys.argv:
    index = sys.argv.index("--Run_experiment")
    if sys.argv[index+1] != "True" and sys.argv[index+1] != "False":
        raise ValueError("Run parameter must be a boolean")
    elif sys.argv[index+1] != "True":
        run_experiment = True

# Check if there is a Time parameter with a positive integer
if "--Time" in sys.argv:
    index = sys.argv.index("--Iterations")
    if type(sys.argv[index+1]) != type(1):
        raise ValueError("Time parameter must be an integer")
    time_end = int(sys.argv[index+1])
    if time_end < 1:
        raise ValueError("Time parameter must be a positive integer")

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
    if run_experiment:
        for algorithm in algorithms:
            start = time.time()
            n_runs = 0
            
            while time.time() - start < time_end:
                print(f"run: {n_runs}")
                one_hour_run(algorithm)
                n_runs += 1
    if type_plot == "one_hour_run":
        visualise_algorithm_data(algorithms, type="one_hour_run")
    elif type_plot == "time":
        visualise_algorithm_data(algorithms, type="time")
    