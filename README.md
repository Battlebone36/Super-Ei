# Protein Pow(d)er
Proteins are long chains of amino acids that play a crucial role in various processes in the human body. In this project, we work with three types of amino acids: hydrophobic amino acids (H), polar amino acids (P), and cysteine amino acids (C).

The stability of a protein is determined by certain interactions between the amino acids:
- Hydrophobic amino acids (H) form a stabilizing bond when they are adjacent to each other.
- Polar amino acids (P) avoid bonds and do not contribute to stability.
- Cysteine amino acids (C) form stabilizing bonds with other cysteines or with hydrophobic amino acids.

The goal is to fold the protein in such a way that the total stability, expressed in the lowest possible score, is maximized.

# Algorithms
To fold the protein in the best possible way we have made multiple algorithms:
- Random
- Depth First
- Greedy
- Hill Climb
- Mountain Climb
- Simulated Annealing
- Genetic

# Installation and setup
To do everything we did: make proteins, fold them with algorithms, visualise them and gather data, the repository needs to be installed and set up. In this paragraph we will walk you through it.
## Requirements
There are some libraries and packages required to run our code like numpy and pandas for calculation and matplotlib and seaborn for visualising.
### Library
The library we used was python.
To install it you need to run the following line, enabling the development of graphical user interfaces (GUIs) in Python on your system. 
```
sudo apt-get install python3-tk
```
### Packages
The requirements.txt file contains all necessary packages to run the code successfully. These can be easily installed via pip using the following instruction:
```
pip install -r requirements.txt
```


## Repository structure
The following list describes the main directories and files in the repository, and where to find them:
* **/code**: Contains all the code for this project.
    * **/code/algorithms**: Contains the code for the algorithms.
    * **/code/classes**: Contains the needed classes for this case.
    * **/code/data**: Contains the data of the experiments in csv files.
    * **/code/visualisation**: contains the code for the visualisation.

* **/static**: Contains static files like pictures and milestone files.
    * **/static/old pictures** Contains old pictures that are not used anymore.
    * **/static/presentation**: A folder for everything that is in the presentation.
        * **/static/presentation/data** Contains the data.
        * **/static/presentation/hillclimb** Contains the pictures of hillclimb iterations.
        * **/static/presentation/protein explenation** Contains pictures of how proteins work.

# Assumptions
Assumptions of specific parameters in the different algorithms.
* Hill climb algorithm
    - Bias for folding at the end of the protein
        The algorithm loops through the protein and stores the folds that are bennificial but if two folds have the best effect, the latter is stored. This means that the algorithm has a bias for folding at the end of the protein
* Simulated annealing
    - Exponential formula
        We use an exponential formula for the cooling, because of previous research in simulated annealing algorithms. Specifically about Rosetta, a protein folding simulation (Wenlong et al., 2006) and a paper where exponential functions were mentioned to solve bio-informatic problems (Kirkpatrick et al., 1983). 
    - Initial-, end temperature
        Initial temperature, end temperature, cooling rate rate and trials per temperature have been found by a looking at multiple sets of values for the parameters.
* Genetic algorithm
    - Mutation probability
    - Mutation per protein

# Usage of the program
To run every file in the root you need to type in the file name and then the Algorithms that you want to use:
`python3 <filename> <Algorithm1> [<Algorithm2> ...]`.\
An example of this is
```
python3 fold_by_algorithm.py HillClimb SimulatedAnnealing
```


## Fold Protein
The file fold_protein_by_algorithm.py can fold a protein according to a given algorithm and will visualise it. The file can be used the following way.\
Usage: `python3 fold_by_algorithm.py <Algorithm1> [<Algorithm2> ...] [--Protein_sequence <str=HHPCHHHPCHPHHCHPH> --Print_Possible_Algorithms <bool=False>]`\
New Arguments:
- `--Protein_sequence <str=HHPCHHHPCHPHHCHPH>`: A sequence of a protein that must consist of P, H and C to fold. 
- `--Print_Possible_Algorithms <bool=False>`: A boolean flag that will print all possible algorithms.
    options:
    - Random
    - Greedy
    - DepthFirst
    - HillClimb
    - MountainClimb
    - SimulatedAnnealing
    - Genetic
An example to run this is:
```
python3 fold_by_algorithm.py HillClimb SimulatedAnnealing --Protein_sequence HHHHPPPPHPPPC
```

## Experiments
### time_experiment
An experiment has been set up to find the best possible solution in a set amount of time. This experiment will run an algorithm for 1 hour and display the results that have been found. The reason for only running 1 hour is that faster algorithms will have a slight edge because they can run more times but the slower algorithms should compensate for this through their better performance.\
How to run the time experiment:\
Usage: `python3 time_experiment.py <Algorithm1> [<Algorithm2> ...] [--Print_Possible_Algorithms <bool=False> --Type_plot <str=one_hour_run>]`\
New Arguments:
- `Run_experiment <bool=False>`: Runs the experiment and stores the data in a file.
- `Print_Possible_Algorithms <bool=False>`: prints all Algorithm options to visualise.
    options:
    - Random
    - Greedy
    - HillClimb
    - MountainClimb
    - SimulatedAnnealing
    - Genetic
- `Type_plot <str=one_hour_run>`: The type of plot that needs to be visualised.
    options:
    - one_hour_run (default)
    - time
An example to run this file is:
```
python3 time_experiment.py HillClimb SimulatedAnnealing --Run_experiment False --Times 500
```
To gather data you need to switch `Run_experiment` to True

### run_many_times_experiment
Another experiment that has been set up is one that runs all the algorithms a certain amount of times. This will result in a distribution of scores which show the occurence of higher and lower scores. Some algorithms might take longer to run than others.
How to run the run experiment:\
Usage: `Usage: python3 run_many_times_experiment.py <Algorithm1> [<Algorithm2> ...] [--Run_experiment <bool=False> --Times <int=100> --Print_Possible_Algorithms <bool=False>]`
New Arguments:
- `Run_experiment <bool=False>`: Runs experiment and stores the data in a file
- `Times <int=100>`: The amount of times that every Algorithms must be run
- `Print_Possible_Algorithms <bool=False>`: prints all Algorithm options to visualise
    options:
    - Random
    - Greedy
    - HillClimb
    - MountainClimb
    - SimulatedAnnealing
    - Genetic
An example to run this is
```
python3 run_many_times_experiment.py HillClimb SimulatedAnnealing --Run_experiment False --Times 500
```
To gather data you need to switch `Run_experiment` to True

### iteration_experiment
Another experiment that has been set up is one that runs all the algorithms with a certain amount of iterations. An iteration is a single state that is being checked. A state is a configuration of the protein. This will result in a graph which shows the amount of steps next to the stability. Allowing for the insight of a how the solutions becomes better over different states.\
How to run the iteration experiment:\
Usage: `python3 iteration_experiment.py <Algorithm1> [<Algorithm2> ...] [--Run_experiment <bool=False> --Times <int=5000> --Print_Possible_Algorithms <bool=False>]`\
New Arguments:
- `Run_experiment <bool=False>`: Runs experiment and stores data in a file
- `Times <int=100>`: The amount of times that every Algorithms must be run
- `Print_Possible_Algorithms <bool=False>`: prints all Algorithm options to visualise
    options:
    - HillClimb
    - MountainClimb
    - SimulatedAnnealing
    - Genetic

An example to run this is
```
python3 run_many_times_experiment.py HillClimb SimulatedAnnealing --Run_experiment False --Times 500
```
To gather data you need to switch `Run_experiment` to True

# Authors
* Sydney Celie
* Bryan Panken
* Ronan Philips

# Bibliography
Li, W., Wang, T., Li, E., Baker, D., Jin, L., Ge, S., ... & Zhang, Y. (2006, April). Parallelization and performance characterization of protein 3D structure prediction of Rosetta. In Proceedings 20th IEEE International Parallel & Distributed Processing Symposium (pp. 8-pp). IEEE. DOI: 10.1109/IPDPS.2006.1639296

Kirkpatrick, S., Gelatt Jr, C. D., & Vecchi, M. P. (1983). Optimization by simulated annealing. science, 220(4598), 671-680. DOI:10.1126/science.220.4598.671