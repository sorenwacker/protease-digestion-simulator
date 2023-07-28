from .Protease import Protease



class Trypsin(Protease):
    def __init__(self):
        super().__init__(name="Trypsin", cleavage_residues=['K', 'R'], no_cleavage_after=['P'], cleavage_position='C')


class Chymotrypsin(Protease):
    def __init__(self):
        super().__init__(name="Chymotrypsin", cleavage_residues=['F', 'Y', 'W'], no_cleavage_after=['P'], cleavage_position='C')


class Pepsin(Protease):
    def __init__(self):
        super().__init__(name="Pepsin", cleavage_residues=['F', 'L'], no_cleavage_after=[], cleavage_position='C')


class Elastase(Protease):
    def __init__(self):
        super().__init__(name="Elastase", cleavage_residues=['A', 'V', 'L'], no_cleavage_after=['P'], cleavage_position='C')


class Thrombin(Protease):
    def __init__(self):
        super().__init__(name="Thrombin", cleavage_residues=['R'], no_cleavage_after=['P'], cleavage_position='C')   
                

class Plasmepsin(Protease):
    def __init__(self):
        super().__init__(name="Plasmepsin", cleavage_residues=['F', 'Y', 'W', 'M', 'L'], no_cleavage_after=[], cleavage_position='C')


class Falcipain2(Protease):
    def __init__(self):
        super().__init__(name="Falcipain2", cleavage_residues=['K', 'R', 'F', 'L'], no_cleavage_after=[], cleavage_position='C')


class Falcipain3(Protease):
    def __init__(self):
        super().__init__(name="Falcipain3", cleavage_residues=['K', 'R', 'F', 'L'], no_cleavage_after=[], cleavage_position='C')


class PfSUB1(Protease):
    def __init__(self):
        super().__init__(name="PfSUB1", cleavage_residues=['R', 'K', 'L'], no_cleavage_after=[], cleavage_position='C')    


# Available proteases dictionary
available_proteases = {
    "Trypsin": Trypsin,
    "Chymotrypsin": Chymotrypsin,
    "Pepsin": Pepsin,
    "Elastase": Elastase,
    "Thrombin": Thrombin,
    "Plasmepsin": Plasmepsin,
    "Falcipain2": Falcipain2,
    "Falcipain3": Falcipain3,
    "PfSUB1": PfSUB1,
}
