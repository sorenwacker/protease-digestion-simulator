import csv

from .PeptideNode import PeptideNode



def generate_peptide_tree(node, proteases, depth=0, max_depth=None, min_length=0):
    if depth >= max_depth:
        return

    for protease in proteases:
        cleaved_peptides = protease.cleave(node.peptide)
        for peptide in cleaved_peptides:
            # create a child node only when the peptide is different from the parent peptide
            if peptide != node.peptide:
                if len(peptide) > min_length:
                    child_node = PeptideNode(peptide, parent=node)
                    node.add_child(child_node)
                    generate_peptide_tree(child_node, proteases, depth + 1, max_depth)
                    
def find_peptide_positions(sequence, peptide):
    positions = []
    index = 0
    while index < len(sequence):
        index = sequence.find(peptide, index)
        if index == -1:
            break
        positions.append(index)
        index += 1
    return positions
    
def print_aligned_peptide(peptide, position, sequence):
    aligned_peptide = [' '] * len(sequence)
    aligned_peptide[position:position + len(peptide)] = peptide
    return ''.join(aligned_peptide)


def print_tree(node, sequence, printed_peptides=None, start_index=0):
    if printed_peptides is None:
        printed_peptides = set()

    if node.parent is not None:
        peptide_position = (node.peptide, start_index)
        if peptide_position not in printed_peptides:
            aligned_peptide = ' ' * start_index + node.peptide
            print(aligned_peptide)
            printed_peptides.add(peptide_position)

    for child in node.children:
        for match in re.finditer(re.escape(child.peptide), sequence[start_index:]):
            child_start_index = match.start() + start_index
            print_tree(child, sequence, printed_peptides, child_start_index + 1)

            
def extract_peptide_sequences(node, root=True):
    peptide_sequences = set()

    if not root:
        peptide_sequences.add(node.peptide)
    
    for child in node.children:
        child_peptide_sequences = extract_peptide_sequences(child, root=False)
        peptide_sequences.update(child_peptide_sequences)
        
    return peptide_sequences


def calculate_possible_cleavage_sites(protease, sequence):
    """
    Calculates the possible number of cleavage sites in the given sequence for the given protease.
    
    Parameters:
    protease (Protease): The protease to check for cleavage sites.
    sequence (str): The peptide sequence to check for cleavage sites.
    
    Returns:
    int: The number of possible cleavage sites in the sequence.
    """
    possible_sites = 0
    for i in range(len(sequence) - 1):
        if sequence[i] in protease.cleavage_residues and sequence[i + 1] not in protease.no_cleavage_after:
            possible_sites += 1
    return possible_sites


def analyze_cleavage_sites(peptide_sequences):
    """
    Analyzes the given peptide sequences for cleavage sites.
    
    Parameters:
    peptide_sequences (list of str): The list of peptide sequences to analyze.
    
    Returns:
    list of str: The list of cleavage sites found in the peptide sequences.
    """
    cleavage_sites = []
    for sequence in peptide_sequences:
        for i in range(len(sequence) - 1):
            cleavage_site = sequence[i:i+2]
            cleavage_sites.append(cleavage_site)
    return cleavage_sites


def identify_known_proteases(peptide_sequences, proteases, original_sequence):
    """
    Identifies known proteases based on the given peptide sequences and original sequence.
    
    Parameters:
    peptide_sequences (list of str): The list of peptide sequences.
    proteases (list of Protease): The list of known proteases.
    original_sequence (str): The original peptide sequence.
    
    Returns:
    list of tuple: Each tuple contains the name of the protease, the number of matched sites, and the probability of match.
    """
    identified_proteases = []

    # Remove duplicates from peptide_sequences
    peptide_sequences = list(set(peptide_sequences))

    for protease in proteases:
        matched_sites_count = 0
        possible_sites_count = calculate_possible_cleavage_sites(protease, original_sequence)
        cleaved_peptides = protease.cleave(original_sequence)
        
        for peptide_sequence in peptide_sequences:
            if peptide_sequence in cleaved_peptides:
                matched_sites_count += 1

        probability = matched_sites_count / possible_sites_count if possible_sites_count != 0 else 0

        identified_proteases.append((protease.name, matched_sites_count, probability))

    return identified_proteases


def print_tree(node, sequence, min_length, printed_peptides=None, start_index=0):
    if printed_peptides is None:
        printed_peptides = set()

    if node.parent is None:
        print(sequence)

    if node.parent is not None:
        peptide_position = (node.peptide, node.parent.peptide)
        if peptide_position not in printed_peptides:
            aligned_peptide = print_aligned_peptide(node.peptide, start_index, sequence)
            if len(node.peptide) < min_length:
                aligned_peptide = '\033[31m' + aligned_peptide + '\033[0m'  # Print in red color
            print(aligned_peptide)
            printed_peptides.add(peptide_position)

    for child in node.children:
        for match in re.finditer(re.escape(child.peptide), sequence[start_index:]):
            child_start_index = match.start() + start_index
            print_tree(child, sequence, min_length, printed_peptides, child_start_index + 1)


def print_aligned_peptide(peptide, start_index, sequence):
    aligned_peptide = [' '] * len(sequence)
    aligned_peptide[start_index - 1:start_index - 1 + len(peptide)] = peptide
    return ''.join(aligned_peptide)    


def write_proteins_to_mascot_csv_file(protein_sequences, acc='alpha', file_name="output.csv"):
    header = ["prot_hit_num", "prot_acc", "prot_mass", "pep_query", "pep_rank", 
              "pep_isbold", "pep_isunique", "pep_exp_mz", "pep_exp_mr", "pep_exp_z",
              "pep_calc_mr", "pep_delta", "pep_miss", "pep_score", "pep_expect",
              "pep_res_before", "pep_seq", "pep_res_after"]

    # default values for other parameters
    pep_query = -1
    pep_rank = -1
    pep_isbold = -1
    pep_isunique = 0
    pep_exp_mz = -1
    pep_exp_mr = -1
    pep_exp_z = -1
    pep_calc_mr = -1
    pep_delta = -1
    pep_miss = -1
    pep_score = -1
    pep_expect = -1
    prot_mass = 0
    pep_res_before = "B"
    pep_res_after = "B"

    with open(file_name, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([])
        writer.writerow(["Protein hits","--------------------------------------------------------"])
        writer.writerow([])        
        writer.writerow(header)
        
        for i, sequence in enumerate(protein_sequences, start=1):
            row = [i, acc, prot_mass, pep_query, pep_rank, pep_isbold, pep_isunique,
                   pep_exp_mz, pep_exp_mr, pep_exp_z, pep_calc_mr, pep_delta, pep_miss,
                   pep_score, pep_expect, pep_res_before, sequence, pep_res_after]
            writer.writerow(row)
        writer.writerow([])            
        writer.writerow(["Peptide matches not assigned to protein hits","--------------------------------------------------------"])
        writer.writerow([]) 