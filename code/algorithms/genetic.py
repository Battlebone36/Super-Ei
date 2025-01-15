from code.classes.protein import Protein
import copy
from code.visualisation.visualisation import *
import random
random.seed(0)
# import numpy as np

def genetic(protein: Protein) -> None:
    proteins_dict: dict[int, tuple[Protein, list[int]]] = {}
    for i in range(2):
        folds = [random.randint(0, 6) for i in range(len(protein.sequence) - 2)]
        proteins_dict[i] = (Protein(protein.sequence), folds)
        for j, DNA in enumerate(folds):
            ans = proteins_dict[i][0].fold_by_DNA(DNA=DNA, index=j + 1)
            possible_folds = {0, 1, 2, 3, 4, 5, 6}
            while (not proteins_dict[i][0].fold_by_DNA(DNA=DNA, index=j + 1) and possible_folds):
            # if ans is False:
                         
                impossible_folds = set()
                impossible_folds.add(DNA)
                possible_folds = possible_folds - impossible_folds
                DNA = random.choice(list(possible_folds))

                print("oh nee", impossible_folds, DNA)



    # print(proteins_dict)
    # visualise_protein(proteins_dict[0])
    # visualise_protein(proteins_dict[1])
    # visualise_protein(proteins_dict[2])
    


test = Protein("HHPPHCPHHC")
genetic(test)