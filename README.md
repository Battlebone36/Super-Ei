# Protein Pow(d)er
Proteins are long chains of amino acids that play a crucial role in various processes in the human body. In this project, we work with three types of amino acids: hydrophobic amino acids (H), polar amino acids (P), and cysteine amino acids (C).

The stability of a protein is determined by certain interactions between the amino acids:
- Hydrophobic amino acids (H) form a stabilizing bond when they are adjacent to each other.
- Polar amino acids (P) avoid bonds and do not contribute to stability.
- Cysteine amino acids (C) form stabilizing bonds with other cysteines or with hydrophobic amino acids.

The goal is to fold the protein in such a way that the total stability, expressed in the lowest possible score, is maximized.

# Installation and setup
## Requirements
### Library
The following command installs the necessary components for Tkinter on Python 3, enabling the development of graphical user interfaces (GUIs) in Python on your system. 
```
sudo apt-get install python3-tk
```
### Packages
The requirements.txt file contains all necessary packages to run the code successfully. These can be easily installed via pip using the following instruction:
```
pip install -r requirements.txt
```

## Usage of the program
[Needs to be changed to the actual run commands]: #
To run the program, use the following command:
```
python3 -m code.visualisation.visualisation
```

## Repository structure
The following list describes the main directories and files in the repository, and where to find them:
* **/code**: contains all the code for this project.
    * **/code/algorithms**: contains the code for the algorithms.
    * **/code/classes**: contains the needed classes for this case.
    * **/code/data**: contains the data of the experiments in csv files
    * **/code/visualisation**: contains the code for the visualisation.
* **/static**: ???

# Assumptions
Assumptions of specific parameters in the different algorithms. --> uitleggen hoe en wat (waarom deze waarden)
EXPLAIN THE CHOICES AND ASSUMPTIONS!!!
* Hill climb algorithm
    - Bias for folding late on in the protein (???)
* Simulated annealing
    - Exponential formula
        We use an exponential formula for the cooling, because of previous research in simulated annealing algorithms. Specifically about Rosetta, a protein folding simulation (Wenlong et al., 2006) and a paper where exponential functions were mentioned to solve bio-informatic problems (Kirkpatrick et al., 1983). 
    - Initial- and end temperature
* Genetic algorithm
    - Mutation probability
    - Mutation per protein instead of per amino acid

# 

# Experiment
## time_experiment
An experiment has been set up to find the best possible solution in a set amount of time. This experiment will run an algorithm for 1 hour and display the results that have been found. The reason for only running 1 hour is that faster algorithms will have a slight edge because they can run more times but the slower algorithms should compensate for this through their better performance.\
How to run the time experiment:\
Usage: python3 time_experiment.py \<Algorithm1> [\<Algorithm2> ...] [--Run_experiment \<bool=False> --Times \<int=100> --Print_Possible_Algorithms \<bool=False>]\
Arguments:
- \<Algorithm1> [\<Algorithm2> ...]: one or more algorithms
- flaggs:
    - Print_Possible_Algorithms \<bool=False>: prints all Algorithm options to visualise.\
        options:
        - Random
        - Greedy
        - HillClimb
        - MountainClimb
        - SimulatedAnnealing
        - Genetic
    - Type_plot \<str=one_hour_run>: The type of plot that needs to be visualised.\
        options:
        - one_hour_run
        - time

## run_many_times_experiment
Another experiment that has been set up is one that runs all the algorithms a certain amount of times. This will result in a distribution of scores which show the occurence of higher and lower scores. Some algorithms might take longer to run than others.
How to run the run experiment:\
Usage: Usage: python3 run_many_times_experiment.py \<Algorithm1> [\<Algorithm2> ...] [--Run_experiment \<bool=False> --Times \<int=100> --Print_Possible_Algorithms \<bool=False>]\
Arguments:
- \<Algorithm1> [\<Algorithm2> ...]: one or more algorithms
    - flaggs:
        - Run_experiment \<bool=False>: Runs experiment
        - Times \<int=100>: The amount of times that every Algorithms must be run
        - Print_Possible_Algorithms \<bool=False>: prints all Algorithm options to visualise
            options:
            - Random
            - Greedy
            - HillClimb
            - MountainClimb
            - SimulatedAnnealing
            - Genetic

## iteration_experiment
Another experiment that has been set up is one that runs all the algorithms with a certain amount of iterations. An iteration is a single state that is being checked. A state is a configuration of the protein. This will result in a graph which shows the amount of steps next to the stability. Allowing for the insight of a how the solutions becomes better over different states.\
How to run the iteration experiment:\
Usage: python3 iteration_experiment.py \<Algorithm1> [\<Algorithm2> ...] [--Run_experiment \<bool=False> --Times \<int=5000> --Print_Possible_Algorithms \<bool=False>]\
Arguments:
Arguments:
- \<Algorithm1> [\<Algorithm2> ...]: one or more algorithms
    - flaggs:
        - Run_experiment \<bool=False>: Runs experiment
        - Times \<int=100>: The amount of times that every Algorithms must be run
        - Print_Possible_Algorithms \<bool=False>: prints all Algorithm options to visualise
            options:
            - HillClimb
            - MountainClimb
            - SimulatedAnnealing
            - Genetic

# Authors
* Sydney Celie
* Bryan Panken
* Ronan Philips

# Bibliography
Li, W., Wang, T., Li, E., Baker, D., Jin, L., Ge, S., ... & Zhang, Y. (2006, April). Parallelization and performance characterization of protein 3D structure prediction of Rosetta. In Proceedings 20th IEEE International Parallel & Distributed Processing Symposium (pp. 8-pp). IEEE. DOI: 10.1109/IPDPS.2006.1639296

Kirkpatrick, S., Gelatt Jr, C. D., & Vecchi, M. P. (1983). Optimization by simulated annealing. science, 220(4598), 671-680. DOI:10.1126/science.220.4598.671