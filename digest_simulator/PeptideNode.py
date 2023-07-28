class PeptideNode:
    def __init__(self, peptide, parent=None):
      """
      This class represents a node in a peptide tree.

      Parameters
      ----------
      peptide : str
          The peptide sequence of the node.
      parent : PeptideNode, optional 
          The parent node of the node. The default is None.
      """
      self.peptide = peptide
      self.parent = parent
      self.children = []

    def add_child(self, child):
        self.children.append(child)

    def __str__(self):
        return self.peptide
