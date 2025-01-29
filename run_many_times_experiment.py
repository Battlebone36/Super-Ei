from code.data.data_store import *
import sys
if len(sys.argv) < 2:
    print("Usage: python3 run_many_times_experiment.py <Algorithm1> [<Algorithm2> ...]")
    print("       [--Run_experiment <bool=False> --Iterations <int=100> --Print_Possible_Algorithms <bool=False>]")
    print("Example: python3 run_many_times_experiment.py HillClimb MountainClimb SimulatedAnnealing Genetic --Run_experiment True --Iterations 300")
    sys.exit(1)

if "--Print_Possible_Algorithms" in sys.argv:
    print("Possible algorithms are: Random, Greedy, HillClimb, MountainClimb, SimulatedAnnealing, Genetic")
    sys.exit(0)

run_experiment = False
times = 100
possible_algorithms = ["random", "greedy", "hillclimb", "mountainclimb", "simulatedannealing", "genetic"]
limit_to_algorithms = 1

# Check if there is a Run parameter with a boolean
if "--Run_experiment" in sys.argv:
    index = sys.argv.index("--Run_experiment")
    if sys.argv[index+1] != "True" and sys.argv[index+1] != "False":
        raise ValueError("Run parameter must be a boolean")
    elif sys.argv[index+1] != "True":
        run_experiment = True

# Check if there is an Iterations parameter with a positive integer
if "--Iterations" in sys.argv:
    index = sys.argv.index("--Iterations")
    if type(sys.argv[index+1]) != type(1):
        raise ValueError("Iterations parameter must be an integer")
    times = int(sys.argv[index+1])
    if times < 1:
        raise ValueError("Iterations parameter must be a positive integer")

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

    # algorithms = [Random,
    #               HillClimb,
    #               MountainClimb,
    #               SimulatedAnnealing,
    #               Genetic]
    
    if run_experiment:
        for algorithm in algorithms:
            for i in range(100):
                store_data(algorithm)

    visualise_algorithm_data(algorithms)