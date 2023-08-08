from itertools import combinations
import pandas as pd

from .DigestionSimulator import DigestionSimulator
from itertools import combinations
import pandas as pd

class ProteasePredictor:
    def __init__(self, original_sequence, proteases, lambda_penalty=0.5, min_peptide_length=3):
        """
        Initialize the ProteasePredictor with the given sequence and proteases.
        
        Parameters:
        - original_sequence (str): Original protein sequence that was digested.
        - proteases (list): List of protease objects to be considered in the prediction.
        - lambda_penalty (float, optional): Penalty weight for unmatched predicted peptides. Default is 0.5.
        - min_peptide_length (int, optional): Minimum peptide length to consider in the prediction. Default is 3.
        """
        self.original_sequence = original_sequence
        self.proteases = proteases
        self.lambda_penalty = lambda_penalty
        self.min_peptide_length = min_peptide_length

    def _simulated_cleave(self, protease_combination):
        """Helper method to simulate cleavage of sequence with a combination of proteases."""
        simulator = DigestionSimulator(self.original_sequence, protease_combination, self.min_peptide_length)
        return set(simulator.extract_unique_peptide_sequences())

    def predict(self, peptide_sequences):
        """
        Predicts the potential proteases responsible for generating observed peptide sequences.
        
        Parameters:
        - peptide_sequences (list): List of peptide sequences observed after digestion.
        
        Returns:
        - DataFrame: A pandas DataFrame sorted by score. Each row contains the combination of proteases,
          the number of matched peptides, and the overall score.
        """
        identified_proteases = []
        peptide_sequences_set = set(peptide_sequences)
        total_peptide_count = len(peptide_sequences_set)

        for i in range(1, len(self.proteases) + 1):  
            for protease_combination in combinations(self.proteases, i):
                names = '+'.join([p.name for p in protease_combination])
                cleaved_peptides_set = self._simulated_cleave(protease_combination)
                
                matched_peptides_count = len(cleaved_peptides_set.intersection(peptide_sequences_set))
                unmatched_predicted_count = len(cleaved_peptides_set) - matched_peptides_count
                probability = matched_peptides_count / total_peptide_count if total_peptide_count else 0
                
                penalty = self.lambda_penalty * unmatched_predicted_count / len(cleaved_peptides_set)
                score = probability - penalty
                
                identified_proteases.append((names, matched_peptides_count, score))

        # Convert to a DataFrame
        df = pd.DataFrame(identified_proteases, columns=['Protease', 'Matched_Peptides', 'Score'])
        
        # Sort by Score in descending order
        df = df.sort_values(by='Score', ascending=False).reset_index(drop=True)
        
        return df
