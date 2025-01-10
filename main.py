from code.classes.protein import Protein
from code.visualisation.visualisation import Visualise
import csv

if __name__ == "__main__":
    protein_vis = Protein("HHCPPPPH", "manual")
    Visualise.visualise_protein(protein_vis)
    protein_vis.fold((2,0), "right")
    Visualise.visualise_protein(protein_vis)

    with open ('output.csv','w',newline = '') as csvfile:
        my_writer = csv.writer(csvfile, delimiter = ',')
        my_writer.writerows(protein_vis.output())
