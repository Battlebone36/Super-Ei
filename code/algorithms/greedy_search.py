from code.classes.protein import Protein
from code.visualisation.visualisation import *
import random
import copy


def greedy_search_sequence(protein: Protein) -> Protein:
    old_score = 0
    max_score = 0
    old_protein = copy.deepcopy(protein)
    old_protein.fold(*random.choice(protein.possible_folds()))
    
    for i in range(30):
        folds = old_protein.possible_folds()
        random.shuffle(folds)
        possible_choices = []
        for fold in folds:
            old_protein.fold(*fold)
            stab = old_protein.stability()
            # possible_choices.append((stab, fold))
            if max_score > stab:
                max_score = stab
                new_protein = copy.deepcopy(old_protein)
            elif max_score + 3 > stab:
                possible_choices.append((stab, fold))
            old_protein.fold_revers(*fold)
        
        if old_score == max_score and possible_choices:
            old_protein.fold(*possible_choices[0][1])
            new_protein = copy.deepcopy(old_protein)
        
        old_protein = copy.deepcopy(new_protein)
        old_score = max_score

    return new_protein
