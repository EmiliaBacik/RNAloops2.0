#!/usr/bin/python
import vg
import numpy as np
from Bio.PDB import *
from Bio.PDB import PDBList
from pathlib import Path
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB import MMCIF2Dict

from Commons.TetriaryStructuresTools import StructureSelection2, points_to_vector_converter, generate_fragments, remove_HOH_from_model, standardize_model, generate_fragments_w_pseudo
from Commons.Utilities import is_non_zero_file, pairs_generator


def centroid_calculator(list_of_atoms):
    """
    Function to find center for a set of atoms
    :param list_of_atoms: list of atom coordinates
    :return: center point
    """
    length = len(list_of_atoms)
    sum1 = map(sum, zip(*list_of_atoms))
    return [val/length for val in sum1]

def planar_angle_calculator(first_stem,second_stem):
    """
    Function to calculate planar angle between 2 vectors
    :param first_stem: stem coordinates
    :param second_stem: stem coordinates
    :return: value of planar angle between two vectors
    """
    return np.nan_to_num(vg.angle(np.array([first_stem[0], first_stem[1], first_stem[2]]), np.array([second_stem[0], second_stem[1], second_stem[2]])))

def euler_angle_calculator(first_stem,second_stem):

    """
    Function to calculate Euler angles between 2 vectors
    :param first_stem: stem coordinates
    :param second_stem: stem coordinates
    :return: value of euler angles between two vectors
    """
    angle_x = np.nan_to_num(vg.angle(np.array([first_stem[0], first_stem[1], first_stem[2]]), np.array([second_stem[0], second_stem[1], second_stem[2]]), look=vg.basis.x))
    angle_y = np.nan_to_num(vg.angle(np.array([first_stem[0], first_stem[1], first_stem[2]]), np.array([second_stem[0], second_stem[1], second_stem[2]]), look=vg.basis.y))
    angle_z = np.nan_to_num(vg.angle(np.array([first_stem[0], first_stem[1], first_stem[2]]), np.array([second_stem[0], second_stem[1], second_stem[2]]), look=vg.basis.z))
    return [format(angle_x, '.3f'), format(angle_y, '.3f'), format(angle_z, '.3f')]

def takeLast(elem):
    return elem[-1]

def calculate_euler_angles_pairwise(list_of_stem_pairs, structure_name, structure_chains, save_substructure = True):

    list_of_points = []
    list_of_euler_angles = []
    list_of_planar_angles = []
    list_of_residues = []
    base_list_of_residues=[]
    base_list_of_points = []
    conn_list_of_residues=[]

    pdb_data_folder = Path("PDB_files/")
    structure_name = structure_name.lower()

    try:
        structure = get_structure(pdb_data_folder, structure_name)
    except:
        with open('err.txt', 'a+') as f:
            f.write(structure_name)
        return [], [], '', [], base_list_of_points, False

    model = structure[0]
    chain_test = []
    chain = None

    for chain_to_include in structure_chains:
        cid = chain_to_include[0]
        if cid not in model:
            continue  # pomijamy pseudo-strandy
        if not chain:
            chain = model[cid]
        chain_test.append(model[cid])

    for pair in list_of_stem_pairs[1]:
        residue_1 = get_residue(chain_test, pair[0])
        residue_2 = get_residue(chain_test, pair[1])

        list_of_residues = [*list_of_residues, [[residue_1.full_id,pair[0]], [residue_2.full_id,pair[1]]]]

        atoms_to_include = [*Selection.unfold_entities(residue_1, 'A'), *Selection.unfold_entities(residue_2, 'A')]
        list_of_atoms = [[atom.get_vector()[0],atom.get_vector()[1],atom.get_vector()[2]] for atom in atoms_to_include]

        list_of_points = [*list_of_points, centroid_calculator(list_of_atoms)]

    for base_pair in list_of_stem_pairs[0]:
        base_residue_1 = get_residue(chain_test, base_pair[0])
        base_residue_2 = get_residue(chain_test, base_pair[1])

        base_list_of_residues = [*base_list_of_residues, [base_residue_1.full_id, base_residue_2.full_id]]

        base_atoms_to_include = [*Selection.unfold_entities(base_residue_1, 'A'), *Selection.unfold_entities(base_residue_2, 'A')]     
        base_list_of_atoms = [[atom.get_vector()[0],atom.get_vector()[1],atom.get_vector()[2]] for atom in base_atoms_to_include]

        base_list_of_points = [*base_list_of_points, centroid_calculator(base_list_of_atoms)]

    center_of_junction = centroid_calculator(base_list_of_points)

    for bp in list_of_stem_pairs[3]:
        if bp:
            conn_residue_1 = get_residue(chain_test, bp[0])
            conn_residue_2 = get_residue(chain_test, bp[1])

            conn_list_of_residues = [*conn_list_of_residues, [[conn_residue_1.full_id ,conn_residue_1.resname], [conn_residue_2.full_id],conn_residue_2.resname]]
        else:
            conn_list_of_residues = [*conn_list_of_residues, ['','']]

    list_of_residues_to_save = []
    for pair in list_of_stem_pairs[2]:
        res_1 = get_residue(chain_test, pair[0])
        res_2 = get_residue(chain_test, pair[1])
        list_of_residues_to_save = [*list_of_residues_to_save, [[res_1.full_id,pair[0]], [res_2.full_id,pair[1]]]]

    filename = []
    #substructure save
    if save_substructure:
        all_included_residues = []
        frag = []
        pairs = generate_fragments_w_pseudo(list_of_residues_to_save)

        for pair in pairs:
            for x in range(pair[0][1], pair[1][1]):
                all_included_residues.append(get_residue(chain_test, x, True))
            all_included_residues.append(get_residue(chain_test, pair[1][1], True))
            opening_residue = get_residue(chain_test, pair[0][1], True)
            closing_residue = get_residue(chain_test,  pair[1][1], True)
            filename = [*filename, [str(opening_residue[2])+ '-' +str(opening_residue[3][1]) , str(closing_residue[2])+ '-' +str(closing_residue[3][1])]]

        name_of_file = save_structure(structure, filename, structure_name, all_included_residues)

    else:
        name_of_file = ''

    #calculation of euler and planar angles
    for pair in pairs_generator(range(len(list_of_points))):
        first_stem = points_to_vector_converter(center_of_junction, list_of_points[pair[0]])
        second_stem = points_to_vector_converter(center_of_junction, list_of_points[pair[1]])

        planar_angle = planar_angle_calculator(first_stem,second_stem)
        list_of_planar_angles.append(format(planar_angle, '.3f'))

        list_of_euler_angles.append(euler_angle_calculator(first_stem,second_stem))

    fragments_to_range = []
    for pair in list_of_stem_pairs[3]:
        if pair:
            start_res = get_residue(chain_test, pair[0], True)
            start = str(start_res[2]) +'-'+ str(start_res[3][1])
            end_res = get_residue(chain_test, pair[1], True)
            end = str(end_res[2]) +'-'+ str(end_res[3][1])
            fragments_to_range.append([start,end])
        else:
            fragments_to_range.append(['',''])

    return list_of_euler_angles, list_of_planar_angles, name_of_file, conn_list_of_residues, fragments_to_range, True

def save_structure(structure, list_of_stem_pairs, structure_PDB_ID, list_of_residues):
    """
    Function creates the file with junction 3D model (*.cif format) and returns its name
    :param structure: the whole PDB structure
    :param list_of_stem_pairs: list of all pairs that junction consist of
    :param structure_PDB_ID: PDB ID of whole structure
    :param list_of_residues: list of residues full IDs that junction consist of
    :return: name of the created file
    """
    list_of_fragments = list_of_stem_pairs
    #list_of_fragments = [[item[0], item[1]-1] for item in list_of_fragments]
    stems_location = '_'.join(str(item) for sublist in list_of_fragments for item in sublist)
    name_of_file = structure_PDB_ID + '_' + str(len(list_of_stem_pairs)) + 'WJ' + '_' + stems_location + '.cif'
    if not is_non_zero_file('./output/structures/' + (structure_PDB_ID + ".cif")):
        io=MMCIFIO()
        io.set_structure(structure)

        name_of_file = structure_PDB_ID + '_' + str(len(list_of_stem_pairs)) + 'WJ' + '_' + stems_location + '.cif'
        if len(name_of_file)>255:
            name_of_file = name_of_file[0:251] + '.cif'
        io.save('./output/structures/' + name_of_file, StructureSelection2(list_of_residues))
    return name_of_file

def get_structure(pdb_data_folder, structure_PDB_ID):
    """
    Function to retrieve information about structure of specified molecule
    :param pdb_data_folder: path to folder contaiting all 3D structures
    :param structure_PDB_ID: structures PDB ID
    :return: the structure object
    """
    if (not is_non_zero_file(pdb_data_folder / (structure_PDB_ID + ".pdb"))) and (not is_non_zero_file(pdb_data_folder / (structure_PDB_ID + ".cif"))):
        PDBList().retrieve_pdb_file(structure_PDB_ID, pdir=pdb_data_folder)

    if is_non_zero_file(pdb_data_folder / (structure_PDB_ID + ".pdb")):
        parser_pdb = PDBParser()
        structure = parser_pdb.get_structure(structure_PDB_ID, pdb_data_folder / (structure_PDB_ID + ".pdb"))
    if is_non_zero_file(pdb_data_folder / (structure_PDB_ID + ".cif")):
        try:
            parser_cif = MMCIFParser(QUIET=True)
            structure = parser_cif.get_structure(structure_PDB_ID, pdb_data_folder / (structure_PDB_ID + ".cif"))
        except:
            PDBList().retrieve_pdb_file(structure_PDB_ID, pdir=pdb_data_folder, file_format='pdb')
            parser_pdb = PDBParser()
            structure = parser_pdb.get_structure(structure_PDB_ID, pdb_data_folder / (structure_PDB_ID + ".pdb"))
    return structure

def get_residue(list_of_chains, res_num, get_full_id = False):
    """
    Function to retrieve the resudie object or residue full id
    :param list_of_chains: list of all avaiable chains
    :param res_num: number of residue to retrieve
    :param get_full_id: bolean if True the function returns residue_full_id, else returns residue object
    :return: residue_full_id/residue object
    """
    sum_of_chain_length = len(list_of_chains[0].child_list)
    chain_id = 0
    while res_num>sum_of_chain_length:
        chain_id += 1
        sum_of_chain_length += len(list_of_chains[chain_id].child_list)
    if get_full_id:
        return list_of_chains[chain_id].child_list[res_num-1-(sum_of_chain_length - len(list_of_chains[chain_id].child_list))].full_id 
    else:
        x= res_num - 1 - (sum_of_chain_length - len(list_of_chains[chain_id].child_list))
        return list_of_chains[chain_id].child_list[res_num-1-(sum_of_chain_length - len(list_of_chains[chain_id].child_list))]

