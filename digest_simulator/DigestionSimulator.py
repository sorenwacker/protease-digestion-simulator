import pandas as pd

from .PeptideNode import PeptideNode
from .proteases import available_proteases
from .tools import generate_peptide_tree, draw_tree, extract_peptide_sequences



class DigestionSimulator:
    def __init__(self, sequence, proteases=None, min_peptide_length=3, min_length_color=5, max_depth=100):
        self.sequence = sequence
        self.min_peptide_length = min_peptide_length
        self.min_length_color = min_length_color
        self.max_depth = max_depth
        self.root = PeptideNode(sequence)
        self.unique_peptide_sequences = None
        self.proteases = proteases
        self.generate_peptide_tree()

    def generate_peptide_tree(self):
        generate_peptide_tree(self.root, self.proteases, 0, max_depth=self.max_depth, min_length=self.min_peptide_length)
        
    def draw_tree(self):
        return draw_tree(self.root, self.sequence, min_length=self.min_length_color, start_index=0)

    def print_tree(self, *args, **kwargs):
        print(self.draw_tree(*args, **kwargs))

    def extract_unique_peptide_sequences(self):
        self.unique_peptide_sequences = extract_peptide_sequences(self.root)
        #print('+'*80)
        #print("Extracted unique peptide sequences (excluding root sequence):")
        #for i, _sequence in enumerate(sorted(self.unique_peptide_sequences, key=len)):
        #    print(i, _sequence)
        #print('+'*80)
        return self.unique_peptide_sequences

