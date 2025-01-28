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
    * **/code/visualisation**: contains the code for the visualisation.
* **/data**: contains the generated csv files and scripts for storing and analyzing algorithm performance data.
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

# Experiment
An experiment has been set up to find the best possible solution in a set amount of time. This experiment will run an algorithm for 1 hour and display the results that have been found. The reason for only running 1 hour is that faster algorithms will have a slight edge because they can run more times but the slower algorithms should compensate for this through their better performance.
How to run the time experiment:
- python3 script1.py

# Authors
* Sydney Celie
* Bryan Panken
* Ronan Philips

# Bibliography
Li, W., Wang, T., Li, E., Baker, D., Jin, L., Ge, S., ... & Zhang, Y. (2006, April). Parallelization and performance characterization of protein 3D structure prediction of Rosetta. In Proceedings 20th IEEE International Parallel & Distributed Processing Symposium (pp. 8-pp). IEEE. DOI: 10.1109/IPDPS.2006.1639296

Kirkpatrick, S., Gelatt Jr, C. D., & Vecchi, M. P. (1983). Optimization by simulated annealing. science, 220(4598), 671-680. DOI:10.1126/science.220.4598.671