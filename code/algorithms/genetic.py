from code.classes.protein import Protein
import copy
from code.visualisation.visualisation import *
import random
# random.seed(0)
# import numpy as np

def genetic(protein: Protein) -> None:
    proteins_dict: dict[int, tuple[Protein, list[int]]] = {}

    # Generation 0
    for i in range(2):
        folds = [random.randint(0, 6) for i in range(len(protein.sequence) - 2)]
        proteins_dict[i] = (Protein(protein.sequence), folds)
        for j, DNA in enumerate(folds):
            proteins_dict[i][0].fold_by_DNA(DNA=DNA, index=j + 1)
            
            # Add stability to the already existing dictionary 
            # Split the two best proteins in two -> and add these together -> based on stability (include weighted options)
                # half 1.1 + half 2.2
                # half 2.1 + half 1.2
            # Keep it going until it can't be done anymore
            



    # print(proteins_dict)
    # visualise_protein(proteins_dict[0])
    # visualise_protein(proteins_dict[1])
    # visualise_protein(proteins_dict[2])
    


test = Protein("HHPPHCPHHC")
genetic(test)