#!/usr/bin/python
import glob
import os
from pathlib import Path

from collections import Counter
import Commons.XML_Generator
from Commons.DataStructures import Record, Stem
from Commons.JunctionFinder import JunctionFinder


def multiple_file_custom_mode(path):
    list_of_junction_types = []
    list_of_records = []

    for filename in glob.glob(os.path.join(path, '*.dbn')):
        name_of_the_file = filename.split("/")[-1]
        name_of_the_file = name_of_the_file.replace('.dbn', '')
        file_info = name_of_the_file.split('-')

        current_record = Record(file_info[0], '', [], 0, [])
        chain_count = 1
        while len(file_info[chain_count]) == 1:
            current_record.chains.append(file_info[chain_count])
            chain_count += 1
        current_record.method = file_info[chain_count]

        current_file = open(filename, "r")

        current_file_content = current_file.readlines()
        sequence = current_file_content[1]
        dot_bracket_representation = current_file_content[2]

        types_of_junction, current_record.list_of_junctions, valid = JunctionFinder.find_junctions(dot_bracket_representation, sequence, file_info[0], current_record.chains, current_record.method)
        list_of_records.append(current_record)
        Commons.XML_Generator.xml_generate_single_rec(current_record)
        list_of_junction_types += types_of_junction

    list_of_records[:] = [dot_bracket_representation for dot_bracket_representation in list_of_records if len(dot_bracket_representation.list_of_junctions) > 0]
    Commons.XML_Generator.xml_generate(list_of_records)
    Commons.XML_Generator.newRunRun()
    #draw_graph_2(list_of_records)
    print(Counter(list_of_junction_types))


def multiple_files_from_pdbee(path):
    list_of_junction_types = []
    list_of_records = []
    for filename in glob.glob(os.path.join(path, '*.dbn')):
        print(filename)
        sequence = ''
        dot_bracket_representation = ''
        name_of_the_file = Path(filename)
        file_info = name_of_the_file.name.split('-')
        if not (os.path.isfile('./output/single_records/' + file_info[0] + '.xml')):
            current_record = Record(file_info[0], '', [], [])
            with open(filename,"r") as dot_bracket_file:
                lines = dot_bracket_file.readlines()
                for line_number in range(len(lines)):
                    if lines[line_number].startswith(">strand"):
                        current_record.chains.append([lines[line_number].split(">strand_")[-1].replace('\n',''),len(lines[line_number + 1].rstrip())])
                        sequence += lines[line_number + 1].rstrip()
                        dot_bracket_representation += lines[line_number + 2].rstrip()
                types_of_junction, current_record.list_of_junctions, valid = JunctionFinder.find_junctions(dot_bracket_representation, 
sequence, file_info[0], current_record.chains, current_record.method)
                if not(valid):
                    continue
                #list_of_records.append(current_record)
                Commons.XML_Generator.xml_generate_single_rec(current_record)
                list_of_junction_types += types_of_junction
        else:
            continue

    #list_of_records[:] = [dot_bracket_representation for dot_bracket_representation in list_of_records if len(dot_bracket_representation.segment_length) > 0]
    #Commons.XML_Generator.xml_generate(list_of_records)
    Commons.XML_Generator.newRunRun()
    print(Counter(list_of_junction_types))


def single_files_from_pdbee(filename):
    list_of_junction_types = []
    list_of_records = []

    sequence = ''
    dot_bracket_representation = ''
    name_of_the_file = Path(filename)
    file_info = name_of_the_file.name.split('-')
    if not (os.path.isfile('./output/single_records/' + file_info[0] + '.xml')):
        current_record = Record(file_info[0], '', [], [])
        with open(filename,"r") as dot_bracket_file:
            lines = dot_bracket_file.readlines()
            for line_number in range(len(lines)):
                if lines[line_number].startswith(">strand"):
                    current_record.chains.append([lines[line_number].split(">strand_")[-1].replace('\n',''),len(lines[line_number + 1].rstrip())])
                    sequence += lines[line_number + 1].rstrip()
                    dot_bracket_representation += lines[line_number + 2].rstrip()
            types_of_junction, current_record.list_of_junctions, valid = JunctionFinder.find_junctions(dot_bracket_representation, sequence, file_info[0], current_record.chains, current_record.method)
            if not(valid):
                return
            #list_of_records.append(current_record)
            Commons.XML_Generator.xml_generate_single_rec(current_record)
            list_of_junction_types += types_of_junction
    else:
        return

    #list_of_records[:] = [dot_bracket_representation for dot_bracket_representation in list_of_records if len(dot_bracket_representation.segment_length) > 0]
    #Commons.XML_Generator.xml_generate(list_of_records)
    #Commons.XML_Generator.newRunRun()
    print(Counter(list_of_junction_types))
