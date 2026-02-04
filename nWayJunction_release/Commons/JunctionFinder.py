#!/usr/bin/python
from Commons.AngularCalculator import calculate_euler_angles_pairwise
from Commons.DataStructures import Junction, Stem, Connector
from Commons.SecondaryStructureTools import find_junction, fill_secondary
from Commons.TetriaryStructuresTools import generate_fragments
from Commons.Utilities import pairs_generator, extract
import itertools
from collections import Counter
import pandas as pd
import copy



class JunctionFinder:
    def __init__(self, text, sequence, name_of_structure, chains, method):
        self.text = text
        self.sequence = sequence
        self.name_of_structure = name_of_structure
        self.chains = chains
        self.method = method

    @classmethod
    def extend_list_of_pairs_w_pseudo1(cls, list_of_pairs):
        list_of_pairs_copy = list_of_pairs
        list_of_pairs = list(itertools.chain(*list_of_pairs))
        list_of_pairs.sort()
        list_of_pairs = list_of_pairs[1::2]
        dict = {}
        for el in list_of_pairs:
            dict[el] = el +1
        for k,c in enumerate(list_of_pairs_copy):
            #if k ==0:
                #list_of_pairs_copy[k][0] = list_of_pairs_copy[k][0]-1
            if list_of_pairs_copy[k][0] in dict.keys():
                list_of_pairs_copy[k][0] = dict[list_of_pairs_copy[k][0]]
            if list_of_pairs_copy[k][1] in dict.keys():
                list_of_pairs_copy[k][1] = dict[list_of_pairs_copy[k][1]]
        return [list_of_pairs_copy]


    @staticmethod
    def find_junctions(text, sequence, name_of_structure, chains, method):

        types_of_junction = []
        list_of_junctions = []
        found_junctions,segments_ranges, order_list = find_junction(text)
        #print("Found junctions: ", found_junctions)
        #print("Segments_ranges: ", segments_ranges)
        #print("Order list: ", order_list)
        if found_junctions:
            for k, j in enumerate(found_junctions):
                number_of_stems = j[0]
                position_of_connectors = j[1]
                opt_max_stem_length = [j[2]]
                current_junction = Junction()
                current_junction.list_of_stems = []
                current_junction.list_of_connectors = []
                current_junction.segment_length =[]
                current_junction.sequence = ''
                current_junction.list_of_segment_db = ''
                lt_cp = copy.deepcopy(position_of_connectors)
                list_of_pairs = [lt_cp]
                extended_of_pairs = copy.deepcopy(position_of_connectors)

                if order_list[k]>0:
                    opt_max_stem_length = [[1 for x in opt_max_stem_length[0]]]
                    list_of_pairs = [*list_of_pairs, lt_cp]
                    list_of_pairs = [*list_of_pairs, lt_cp]
                    #list_of_pairs = [*list_of_pairs, *JunctionFinder.extend_list_of_pairs_w_pseudo1(extended_of_pairs)]
                    #list_of_pairs = [*list_of_pairs, *JunctionFinder.extend_list_of_pairs(opt_max_stem_length, number_of_stems, position_of_connectors)]
                    list_of_pairs = [*list_of_pairs, segments_ranges[k]]

                elif opt_max_stem_length:
                    limited = [[3 if x >= 3 else x for x in opt_max_stem_length[0]]]
                    list_of_pairs = [*list_of_pairs, *JunctionFinder.extend_list_of_pairs(limited, number_of_stems, position_of_connectors)]
                    list_of_pairs = [*list_of_pairs, *JunctionFinder.extend_list_of_pairs(opt_max_stem_length, number_of_stems, position_of_connectors)]
                    list_of_pairs = [*list_of_pairs, segments_ranges[k]]

                euler_collection = []
                planar_collection = []
                files_collection = []

                list_of_euler_angles, planar_angle, junction_substructure_file, id_pairs, fragments, valid= calculate_euler_angles_pairwise(list_of_pairs, name_of_structure, chains)
                if not(valid):
                    return [], [], False
                euler_collection.append(list_of_euler_angles)
                planar_collection.append(planar_angle)
                files_collection.append(junction_substructure_file)

                #fill values for junction
                current_junction.type = str(number_of_stems) + '-way junction'
                current_junction.name_of_file = files_collection

                current_junction.segment_length.append(3)

                list_of_fragments = segments_ranges[k]
                #print("Junction", k)
                #print("Connector positions:", position_of_connectors)
                #print("Fragments:", list_of_fragments)


                segment_ranges_ids = id_pairs

                list_db = []
                list_db, srt_db = fill_secondary(list_of_fragments,text)
                for i, pair in enumerate(list_of_fragments):
                    if i == len(list_of_fragments)-1 and order_list[k] == 0:
                        sequence1=(srt_db[pair[0] - opt_max_stem_length[0][i]:pair[1]+opt_max_stem_length[0][0]-1])
                        list_of_segment_db=(sequence[pair[0]- opt_max_stem_length[0][i]:pair[1]+opt_max_stem_length[0][0]-1])
                    elif i == len(list_of_fragments)-1 and order_list[k] > 0:
                        sequence1=(srt_db[pair[0] - opt_max_stem_length[0][i]:pair[1]])
                        list_of_segment_db=(sequence[pair[0]- opt_max_stem_length[0][i]:pair[1]])

                    elif i == 0 and order_list[k] == 0:
                        sequence1=(srt_db[pair[0] - opt_max_stem_length[0][i]:pair[1]+opt_max_stem_length[0][i+1]-1]) + '-'
                        list_of_segment_db=(sequence[pair[0] - opt_max_stem_length[0][i]:pair[1]+opt_max_stem_length[0][i+1]-1]) + '-' 

                    else:
                        sequence1=(srt_db[pair[0] - opt_max_stem_length[0][i]:pair[1]+opt_max_stem_length[0][i+1]-1]) + '-'
                        list_of_segment_db=(sequence[pair[0] - opt_max_stem_length[0][i]:pair[1]+opt_max_stem_length[0][i+1]-1]) + '-' 
                    current_junction.sequence += sequence1
                    current_junction.list_of_segment_db += list_of_segment_db

                list_of_connectors = [JunctionFinder.set_connector_record(i, segments_ranges[k], segment_ranges_ids, text, pair, sequence, euler_collection, planar_collection, fragments) for i,pair in enumerate(list_of_fragments)]

                connector_pairs = []
                for connector in pairs_generator(range(len(list_of_connectors))):
                    connector_pairs.append(connector)
                connector_pairs.insert(0, connector_pairs.pop())
                connector_pairs1 = []
                pos_list = []
                all_pos = []
                #print("position_of_connectors:", position_of_connectors)
                #print("list_of_fragments:", list_of_fragments)

                for c_id, x in enumerate(position_of_connectors):
                    pos1 = []
                    pos2 = []
                    for c1_id, y in enumerate(list_of_fragments):

                        if(y[0]<=x[0]<=y[1]):
                            pos1.append(c1_id)

                        if(y[0]<=x[1]<=y[1]):
                            pos2.append(c1_id)
                    if len(pos1) ==1 or len(pos2)==1:
                        pos_list= (pos1[0], pos2[0])
                    elif len(pos1) >1 or len(pos2)>1:
                        pos_list= (pos1, pos2)
                    all_pos.append(sorted(pos_list))
                for pz, element in enumerate(all_pos):

                    if type(element[0]) is not int:

                        positions_to_alter = []
                        for pz1, el in enumerate(all_pos):
                            if el == element:

                                positions_to_alter.append(pz1)
                        cnt1= 0
                        cnt2 = len(positions_to_alter)
                        for c in positions_to_alter:
                            all_pos[c] = [all_pos[c][0][cnt1], all_pos[c][1][cnt2-1]]
                            cnt1 +=1
                            cnt2 -=1

                connector_pairs1 = all_pos

                current_junction.list_of_stems = [JunctionFinder.set_stem_record(opt_max_stem_length[0], connector_id, sequence, position_of_connectors, connector_pairs1, order_list[k] ) for connector_id in range(number_of_stems)]
                current_junction.list_of_connectors.append(list_of_connectors)
                list_of_junctions.append(current_junction)
                types_of_junction.append(number_of_stems)

        return types_of_junction, list_of_junctions, True

    @classmethod
    def extend_list_of_pairs(cls, max_stem_length, number_of_stems, position_of_connectors):
        extended_list_of_pairs = []

        current_pair = []
        for i in range(number_of_stems):
            if i == 0:
                start = position_of_connectors[i][0]  - max_stem_length[0][i] + 1 if position_of_connectors[i][0] - max_stem_length[0][i] >=0 else position_of_connectors[i][0]
                end =  position_of_connectors[i][1]  + max_stem_length[0][i] - 1 if position_of_connectors[i][0] - max_stem_length[0][i] >=0 else position_of_connectors[i][1]
            else:
                start = position_of_connectors[i][0] + max_stem_length[0][i] - 1
                end =  position_of_connectors[i][1] - max_stem_length[0][i] + 1
            current_pair.append([start, end])
        extended_list_of_pairs.append(current_pair)
        return extended_list_of_pairs

    @classmethod
    def set_stem_record(cls, max_stem_length, connector_id, sequence, position_of_connectors, connector_pairs, k):
        current_stem = Stem()
        current_stem.segment_length = ''

        current_stem.list_of_connectors = []
        current_stem.sequence = []

        if connector_id == 0 and k == 0:
            start = sequence[position_of_connectors[connector_id][0] - max_stem_length[connector_id]: position_of_connectors[connector_id][0]]
            end = sequence[position_of_connectors[connector_id][1] - 1:position_of_connectors[connector_id][1] + max_stem_length[connector_id] - 1]
        else:
            start = sequence[position_of_connectors[connector_id][0] - 1 :position_of_connectors[connector_id][0] + max_stem_length[connector_id] - 1]
            end = sequence[position_of_connectors[connector_id][1] - max_stem_length[connector_id]: position_of_connectors[connector_id][1]]
        stem_sequence =  start + '-' + end
        current_stem.sequence.append(stem_sequence)
        current_stem.segment_length = len(start)
        current_stem.list_of_connectors.append(connector_pairs[connector_id][0])
        current_stem.list_of_connectors.append(connector_pairs[connector_id][1])
        return current_stem

    @classmethod
    def set_connector_record(cls, i, segments_ranges, segment_ranges_ids, text, pair, sequence, euler_collection, planar_collection, fragments):
        current_connector = Connector()
        current_connector.list_of_angles = []
        current_connector.planar_angle = []
        current_connector.list_of_segments_ranges = fragments[i] if segments_ranges[i][0]<=segments_ranges[i][1] else ''
        current_connector.list_of_segments_ranges_id = segment_ranges_ids[i]
        current_connector.lengths_of_segments=(pair[1] - pair[0]-1)
        current_connector.list_of_segment_db=(text[pair[0] - 1:pair[1]])
        current_connector.list_of_segment_seq=(sequence[pair[0] - 1:pair[1]])
        current_connector.list_of_angles.append(extract(i,euler_collection))
        current_connector.planar_angle.append(extract(i, planar_collection))
        current_connector.connector_id=i
        return current_connector

    def max_le(seq, val):
        max_val =  max([v for v in seq])

