from code.classes.protein import Protein
from code.visualisation.visualisation import *
import random


def greedy_search(protein: Protein) -> Protein:
    old_score = 0
    max_score = 0

    for fold in protein.possible_folds():
        protein.fold(*fold)

        a=1
    print(*fold)
        # new_prot = protein.fold(*fold)

    # print(random.choice(protein.possible_folds()))
