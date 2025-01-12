from code.classes.protein import Protein
from code.visualisation.visualisation import *
import csv

if __name__ == "__main__":
    protein_vis = Protein("HHCPPPPH", "manual")
    # print(protein_vis.possible_folds())
    visualise_protein(protein_vis)
    # protein_vis.fold((3, 0, 0), "z_pos")
    # visualise_protein(protein_vis)
    # protein_vis.fold((3, 0, 0), "z_neg")
    # protein_vis.fold((3, 0, 0), "z_neg")
    # visualise_protein(protein_vis)


    with open ('output.csv','w',newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ',')
        my_writer.writerows(protein_vis.output())
