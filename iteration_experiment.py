from code.data.data_store import *
import sys
if len(sys.argv) < 2:
    print("Usage: python3 iteration_experiment.py <Algorithm1> [<Algorithm2> ...]")
    print("       [--Run_experiment <bool=False> --times <int=5000> --Print_Possible_Algorithms <bool=False>]")
    print("Example: python3 iteration_experiment.py HillClimb MountainClimb SimulatedAnnealing Genetic")
    sys.exit(1)

if "--Print_Possible_Algorithms" in sys.argv:
    print("Possible algorithms are: HillClimb, MountainClimb, SimulatedAnnealing, Genetic")
    sys.exit(0)

run_experiment = False
times = 5000
possible_algorithms = ["hillclimb", "mountainclimb", "simulatedannealing", "genetic"]
limit_to_algorithms = 1

# Check if there is a Run_experiment parameter with a boolean
if "--Run_experiment" in sys.argv:
    index = sys.argv.index("--Run_experiment")
    if sys.argv[index+1] != "True" and sys.argv[index+1] != "False":
        raise ValueError("Run_experiment parameter must be a boolean")
    elif sys.argv[index+1] != "True":
        run_experiment = True

# Check if there is an times parameter with a positive integer
if "--times" in sys.argv:
    index = sys.argv.index("--times")
    if type(sys.argv[index+1]) != type(1):
        raise ValueError("times parameter must be an integer")
    times = int(sys.argv[index+1])
    if times < 1:
        raise ValueError("times parameter must be a positive integer")

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
    if algorithms[i] == "hillclimb":
        algorithms[i] = HillClimb
    elif algorithms[i] == "mountainclimb":
        algorithms[i] = MountainClimb
    elif algorithms[i] == "simulatedannealing":
        algorithms[i] = SimulatedAnnealing
    elif algorithms[i] == "genetic":
        algorithms[i] = Genetic

if __name__ == "__main__":
    
    
    print(algorithms)
    
    if run_experiment:
        for algorithm in algorithms:
            run_for_iterations_and_stability(algorithm, times)

    visualise_algorithm_data(algorithms, type="iteration_stability")