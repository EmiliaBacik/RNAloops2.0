import io
import os
from unittest import TestCase
from pathlib import Path
from shutil import rmtree

from Commons.XML_Generator import xml_generate
from Commons.DataStructures import Junction, Record


class TestXml_generator(TestCase):
    def test_xml_generate(self):
        if not os.path.exists(Path('./output')):
            os.makedirs(Path('output'))
        list_of_records = []
        list_of_junctions = []
        current_record = Record('TST', '', [], 0, [])
        current_record.chains = ['R']

        current_junction = Junction()

        number_of_stems = 4
        current_junction.type = str(number_of_stems) + '-way junction'

        raw_list_of_segments_ranges = [[7, 10], [24, 26], [42, 47], [63, 64]]
        current_junction.list_of_segments_ranges = raw_list_of_segments_ranges

        raw_lengths_of_segments = [2, 1, 4, 0]
        current_junction.lengths_of_segments = raw_lengths_of_segments

        raw_list_of_segment_db = ['(..(', ').(', ')....(', '))']
        current_junction.list_of_segment_db = raw_list_of_segment_db

        raw_list_of_segment_seq = ['UuAA', 'UAG', 'UCUAGU', 'AA']
        current_junction.list_of_segment_seq = raw_list_of_segment_seq

        current_junction.list_of_angles = [[168.76346367903653, 154.95465900623543, 147.55895454582856], [19.830967625382193, 23.374574255577755, 25.288146043342312], [172.8539822216024, 161.94480908787452, 147.13728864865644]]
        current_junction.planar_angle = [156.4728220179626, 27.19084937771678, 161.9182427422165, 18.38522279388394]
        current_junction.name_of_file = 'TST_3-way_junction_7_10_24_26_42_47_63_64.pdb'


        list_of_junctions.append(current_junction)
        current_record.list_of_junctions = list_of_junctions
        list_of_records.append(current_record)
        xml_generate(list_of_records)
        tst_path = Path('./output/RESULTS.xml')
        ref_path = Path('./tests/test_files/TEST_RESULTS.xml')
        self.assertListEqual(
            list(io.open(tst_path)),
            list(io.open(ref_path)))
        rmtree(Path("./output"))
