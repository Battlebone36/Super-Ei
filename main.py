from code.classes.protein import Protein
from code.visualisation.visualisation import *
import csv

if __name__ == "__main__":
    protein_vis = Protein("HHCPPPPH", "manual")
    visualise_protein(protein_vis)
    protein_vis.fold((2,0), "right")
    visualise_protein(protein_vis)
    print(f"{protein_vis.output()}")

    with open('output.csv', 'w', newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter=' ')
        my_writer.writerow(protein_vis.output())
