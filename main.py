from code.classes.protein import Protein
from code.visualisation.visualisation import Visualise

protein_vis = Protein("HHCPPPPH", "manual")
Visualise.visualise_protein(protein_vis)
protein_vis.fold((2,0), "right")
Visualise.visualise_protein(protein_vis)

