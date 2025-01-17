from code.classes.protein import Protein
from code.visualisation.visualisation import *

class Algorithm:
    def __init__(self, protein: Protein):
        self.protein = protein

    def run(self):
        self.visualise()

    def visualise(self):
        visualise_protein(self.protein)