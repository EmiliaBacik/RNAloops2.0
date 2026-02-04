import numpy as np
import io
import os
from unittest import TestCase
from pathlib import Path
from Bio.PDB import *
from shutil import rmtree

from Commons.AngularCalculator import planar_angle_calculator, euler_angle_calculator, centroid_calculator, save_structure


class TestAngularCalculator(TestCase):

    def test_euler_angle_calculator_no_change(self):

        first_stem = [0,0,0]
        second_stem = [0,0,0]

        self.assertEqual(euler_angle_calculator(first_stem, second_stem)[0], "{:.3f}".format(0))
        self.assertEqual(euler_angle_calculator(first_stem, second_stem)[1], "{:.3f}".format(0))
        self.assertEqual(euler_angle_calculator(first_stem, second_stem)[2], "{:.3f}".format(0))

    def test_euler_angle_calculator_first_angle(self):

        first_stem = [0,0,1]
        second_stem =[0,1,0]

        self.assertEqual(euler_angle_calculator(first_stem,second_stem)[0], "{:.3f}".format(90))
        self.assertEqual(euler_angle_calculator(first_stem, second_stem)[1], "{:.3f}".format(0))
        self.assertEqual(euler_angle_calculator(first_stem, second_stem)[2], "{:.3f}".format(0))

    def test_euler_angle_calculator_second_angle(self):

        first_stem = [1,0,0]
        second_stem = [0,0,1]

        self.assertEqual(euler_angle_calculator(first_stem, second_stem)[0], "{:.3f}".format(0))
        self.assertEqual(euler_angle_calculator(first_stem,second_stem)[1], "{:.3f}".format(90))
        self.assertEqual(euler_angle_calculator(first_stem, second_stem)[2], "{:.3f}".format(0))

    def test_euler_angle_calculator_second_angle(self):

        first_stem = [1,0,0]
        second_stem = [0,1,0]

        self.assertEqual(euler_angle_calculator(first_stem, second_stem)[0], "{:.3f}".format(0))
        self.assertEqual(euler_angle_calculator(first_stem, second_stem)[1], "{:.3f}".format(0))
        self.assertEqual(euler_angle_calculator(first_stem,second_stem)[2], "{:.3f}".format(90))

    def test_planar_angle_calculator(self):


        first_stem = [0,0,1]
        second_stem =[0,1,0]

        self.assertEqual(planar_angle_calculator(first_stem,second_stem), 90.0)

    def test_planar_angle_calculator_no_change(self):

        first_stem = [6,6,6]
        second_stem =[1,1,1]

        self.assertEqual(planar_angle_calculator(first_stem,second_stem), 0.0)

    def test_centroid_calculator_cube(self):

        #defines points for cube
        set_of_points_cube = [[0,0,0], [0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1]]
        self.assertEqual(centroid_calculator(set_of_points_cube), [0.5, 0.5, 0.5])

    def test_centroid_calculator_triangle(self):
        #triangle on the same plane
        set_of_points_rectangle = [[0,0,0], [0,3,0], [3,0,0]]
        self.assertEqual(centroid_calculator(set_of_points_rectangle), [1, 1, 0])

    def test_save_structure(self):
        if not os.path.exists(Path('./output/structures')):
            os.makedirs(Path('./output/structures'))
        parser_cif = MMCIFParser()
        pdb_data_folder = Path('./tests/test_files')
        structure = parser_cif.get_structure('1b23', pdb_data_folder / '1b23.cif')
        list_of_residues = [[7, 66], [10, 25], [27, 43], [49, 65]]
        structure_name = '1b23'
        list_of_stem_pairs= [[7, 64], [10, 24], [26, 42], [47, 63]]
        save_structure(structure, list_of_stem_pairs, structure_name, list_of_residues)
        tst_path = Path('./output/structures/1b23_4-way_junction_7_10_24_26_42_47_63_64.cif')
        ref_path = Path('./tests/test_files/1b23_4-way_junction_7_10_24_26_42_47_63_64.cif')
        self.assertListEqual(
            list(io.open(tst_path)),
            list(io.open(ref_path)))
        rmtree(Path("./output"))

