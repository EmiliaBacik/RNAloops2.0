#!/usr/bin/python
import re
from Bio.PDB import *
import itertools


def remove_nonRNA_atoms_from_model(model):
    residue_to_remove = []
    chain_to_remove = []
    for chain in model:
        for residue in chain:
            if str(residue.resname).strip() not in ['A', 'C', 'G', 'U' ]:
                    residue_to_remove.append((chain.id, residue.id))
            if len(chain) == 0:
                chain_to_remove.append(chain.id)

    for residue in residue_to_remove:
        model[residue[0]].detach_child(residue[1])

    for chain in chain_to_remove:
        model.detach_child(chain)
    return model

def remove_HOH_from_model(structure):
    list_of_aa = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLU', 'GLN', 'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    residue_to_remove= []
    chain_to_remove = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if (residue.id[0] in ('W', 'H_SO4')) or residue.resname in list_of_aa:
                    residue_to_remove.append((chain.id, residue.id))

            if len(chain) == 0:
                chain_to_remove.append(chain.id)

    for residue in residue_to_remove:
        model[residue[0]].detach_child(residue[1])

    for chain in chain_to_remove:
        model.detach_child(chain)
    return structure

def standardize_model(structure):
    non_standard_residues_A = ['00A', 'A2M', '0DA', '12A', '1AP', '1MA', '26A', '2BU', '2DA', '2MA',
                               '31H', '31M', '3DA', '45A', '5AA', '5FA', '6HB', '6IA', '6MA', '6MZ',
                               '6NW', '7AT', '7DA', '8AN', '8BA', 'A23', 'A2L', 'A2M', '6MA', 'A40',
                               'A38', 'A2M', 'A3A', 'A3P', 'A40', 'A43', 'A44', 'A47', 'A7E', 'A9Z',
                               'ABR', 'ABS', 'AD2', 'ADI', 'ADP', 'AF2', 'AMD', 'AMO', 'AP7', 'AS',
                               'AVC', 'DA',  'E',   'F3N', 'G3A', 'M7A', 'MA6', 'MA7', '1MA', 'MIA',
                               'PPU', 'PU',  'R',   'RBD', 'RIA', 'RMP', 'S4A', 'SDE', 'SMP', 'SRA',
                               'TBN', 'V3L', 'XUA', 'Y', 'MGQ']
    non_standard_residues_G = ['0AD', '0G',  '18M', '1MG', '23G', '2EG', '2JV', '2MG', '2SG', '5CG',
                               '6OG', '7GU', '7MG', '8AA', '8AG', '8FG', '8MG', '8OG', '8OS', 'B8K',
                               'B8W', 'B9B', 'BGM', 'CG1', 'DGP', 'DDG', 'DG',  'DG8', 'DGI', 'DGP',
                               'E6G', 'E7G', 'EFG', 'EQ4', 'F74', 'G1G', 'G2L', 'G2S', 'G31', '6OG',
                               '8MG', 'G36', 'G38', 'G3A', '8OG', 'G46', 'G47', 'G48', 'G49', 'G7M',
                               'GDO', 'GDP', 'GDR', 'GF2', 'GFL', 'GH3', 'GMS', 'GN7', 'GNG', 'GRB',
                               'GS',  'GSR', 'GSS', 'GTP', 'GX1', 'HN0', 'IG',  'KAG', 'KAK', 'LDG',
                               'LGP', 'M2G', 'MG1', 'MGT', 'MGV', 'MHG', 'MRG', 'O2G', 'OGX', 'OMG',
                               'P',   'P7G', 'PGN', 'PGP', 'PPW', 'QUO', 'RDG', 'RFJ', 'S4G', 'S6G',
                               'TGP', 'TPG', 'XPB', 'XUG', 'YYG', 'IGU']
    list_of_nucleotides = ['A', 'C', 'G', 'U', '1MA', '2MG', '5MC', '5MU', '7MG', 'H2U', 'M2G', 'OMC', 'OMG', 'PSU', 'YG', '4SU', 'MIA', 'I', 'DA', 'DT',
                             'DC', 'DU', 'DI', 'DA', 'DG', 'OHX', 'FHU', 'AET']
    non_standard_residues = {'MIA':'A', '1MA':'A', '2MG':'G', '5MC': 'C', '5MU': 'U', '7MG': 'G', 'H2U': 'U', 'M2G':'G', 'OMC':'C', 'OMG':'G', 'PSU':'U', 'YG':'G', '4SU':'U', 'DA': 'A', 'DT': 'T',
                             'DC': 'C', 'DU': 'U', 'DI': 'I', 'DA': 'A', 'DG:':'G', 'FHU': 'U', 'AET':'A'}
    acceptable_atoms = ['C2', 'C4', 'C6', 'C8', 'N1', 'N2', 'N3', 'N4', 'N6', 'N7', 'N9', 'O2', 'O4', 'O6', 'C1\'', 'C2\'', 'C3\'', 'C4\'', 'C5\'', 'O2\'', 'O3\'', 'O4\'', 'O5\'', 'OP1', 'OP2', 'P' ]
    list_of_aa = ['ALA', 'ARG', 'ASN', 'ASP', 'ASX', 'CYS', 'GLU', 'GLN', 'GLX', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']
    residue_to_remove= []
    chain_to_remove = []
    atom_to_remove = []
    #extended = [pair[0] for pair in dictionary_non_standart_residue]
    #list_of_nucleotides = [*list_of_nucleotides, *extended]
    #non_standard_residues = dict(dictionary_non_standart_residue)
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.resname not in list_of_nucleotides:

                    residue_to_remove.append((chain.id, residue.id))
                if residue.resname in non_standard_residues:

                    residue.resname = non_standard_residues[residue.resname]
                    residue_id = list(residue.id)
                    residue_id[0] = ' '
                    residue.id = tuple(residue_id)
                    for atom in residue:
                        if atom.id not in acceptable_atoms:

                            atom_to_remove.append((chain.id, residue.id, atom.id))
                            atom_id = list(atom.full_id[3])
                            atom_id[0] = ' '
                            atom_tuple = tuple(atom_id)
                            atom.full_id = [atom.full_id[0], atom.full_id[1], atom.full_id[2], atom_tuple, atom.full_id[4]]


            if len(chain) == 0:
                chain_to_remove.append(chain.id)
    for atom in atom_to_remove:
        model[atom[0]][atom[1]].detach_child(atom[2])

    for residue in residue_to_remove:
        print(residue)
        model[residue[0]].detach_child(residue[1])

    for chain in chain_to_remove:
        model.detach_child(chain)
    return structure


def points_to_vector_converter(center_of_junction, atom):
    x = atom[0] - center_of_junction[0]
    y = atom[1] - center_of_junction[1]
    z = atom[2] - center_of_junction[2]
    vector = [x, y, z]
    return vector


def generate_fragments(list_of_pairs):
    fragments = []
    for elements in range(1, len(list_of_pairs)):
        if elements == 1:
            fragments.append(sorted([list_of_pairs[0][0], list_of_pairs[elements][0]]))
        elif elements == len(list_of_pairs)-1:
            fragments.append(sorted([list_of_pairs[elements - 1][1], list_of_pairs[elements][0]]))
            fragments.append(sorted([list_of_pairs[elements][1], list_of_pairs[0][1]]))
        else:
            fragments.append(sorted([list_of_pairs[elements - 1][1], list_of_pairs[elements][0]]))
    return fragments

def generate_fragments_w_pseudo(list_of_pairs):
    fragments = []
    short_pairs = list(itertools.chain(*list_of_pairs))
    short_pairs.sort(key = takeLast)
    for elements in range(0, len(list_of_pairs)):
        fragments.append([short_pairs[elements*2], short_pairs[(elements*2)+1]])
    return fragments

def takeLast(elem):
    return elem[-1]

class StructureSelection(Select):
    def __init__(self, list_of_residues):
        self.list_of_residues = list_of_residues

    def accept_residue(self, residue):
        fragment_list = generate_fragments(self.list_of_residues)
        resudue_number_to_include = []
        list_res = ['A', 'U', 'C', 'G']
        for pair in fragment_list:
            for numbers in range(pair[0], pair[1]+1):
                resudue_number_to_include.append(numbers)
        if (residue.id[1] in resudue_number_to_include) and (re.sub(r'\W+', '',residue.get_resname()) in list_res):
            return 1
        else:
            return 0

class StructureSelection2(Select):
    def __init__(self, list_of_residues):
        self.list_of_residues = list_of_residues

    def accept_residue(self, residue):

        if (residue.full_id in self.list_of_residues):
            return 1
        else:
            return 0
