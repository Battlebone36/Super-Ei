from code.classes.protein import Protein
from code.visualisation.visualisation import *
import random
import copy

def greedy_fold(protein: Protein) -> Protein:
    """Folds the protein according to a greedy algorithm which takes the lowest
    score for each possibility."""
    old_protein = copy.deepcopy(protein)
    choice_dict = {}
    lowest_score = 0
    best_direction = ""
    for i in range(len(protein.data)):
        current_protein = copy.deepcopy(old_protein)
        for coord in current_protein.data:
            if current_protein.data[coord] == ("H", i):
                amino = coord
            elif current_protein.data[coord] == ("P", i):
                amino = coord
            elif current_protein.data[coord] == ("C", i):
                amino = coord
            else:
                continue
            
        possible_folds = current_protein.possible_folds_point(amino)
        
        for p_folds in possible_folds:
            current_protein.fold(amino,p_folds)
            new_stability = current_protein.stability()
            choice_dict[p_folds] = new_stability
            # print(choice_dict)
            current_protein = copy.deepcopy(old_protein)
        
        for choices in choice_dict:
            if choice_dict[choices] <= lowest_score:
                lowest_score = choice_dict[choices]
                best_direction = choices


        old_protein.fold(amino, best_direction)
        choice_dict.clear()
    
    return old_protein


    

# protein1 = Protein("HCPHPCPHPCHCHPHPPPHPPPHPPPPHPCPHPPPHPHHHCCHCHCHCHH")
# protein1 = greedy_fold(protein1)
# visualise_protein(protein1)
