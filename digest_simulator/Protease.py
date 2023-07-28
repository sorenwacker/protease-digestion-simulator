class Protease:
    def __init__(self, name, cleavage_residues, no_cleavage_after, cleavage_position):
        """
        This class represents a protease.

        Parameters
        ----------
        name : str
            The name of the protease.
        cleavage_residues : list
            A list of residues that the protease cleaves after.
        no_cleavage_after : list
            A list of residues that the protease does not cleave after.
        cleavage_position : str
            The position of the cleavage. Use 'C' for C-terminal or 'N' for N-terminal.
        """
        self.name = name
        self.cleavage_residues = cleavage_residues
        self.no_cleavage_after = no_cleavage_after
        self.cleavage_position = cleavage_position

    def cleave(self, sequence):
        peptides = []
        start = 0

        if self.cleavage_position == 'C':
            for i in range(len(sequence) - 1):
                if sequence[i] in self.cleavage_residues and sequence[i + 1] not in self.no_cleavage_after:
                    peptides.append(sequence[start:i + 1])
                    start = i + 1
        elif self.cleavage_position == 'N':
            for i in range(1, len(sequence)):
                if sequence[i] in self.cleavage_residues and sequence[i - 1] not in self.no_cleavage_after:
                    peptides.append(sequence[start:i])
                    start = i
        else:
            raise ValueError("Invalid cleavage_position value. Use 'C' for C-terminal or 'N' for N-terminal.")
        
        peptides.append(sequence[start:])
        return peptides
    


