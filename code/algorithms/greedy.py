from code.classes.protein import Protein
from code.visualisation.visualisation import *
import random

def greedy_fold(protein: Protein) -> Protein:
    """folds the protein according to a greedy algorithm which takes the lowest
    score for each possibility"""
    