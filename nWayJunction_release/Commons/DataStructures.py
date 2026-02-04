#!/usr/bin/python


class Record:
    def __init__(self, name, method='', chains=[], list_of_junctions=[]):
        self.name = name
        self.method = method
        self.chains = chains
        self.list_of_junctions = list_of_junctions


class Stem:
    def __init__(self, segment_length=1, list_of_connectors=[]):
        self.segment_length = None
        self.list_of_connectors = list_of_connectors
        self.sequence = None


class Junction:
    def __init__(self):
        self.type = ''
        self.name_of_file = None
        self.list_of_stems = None
        self.list_of_connectors = None
        self.segment_length = None
        self.sequence = None
        self.list_of_segment_db = None

class Connector:
     def __init__(self):
        self.connector_id = None
        self.list_of_segments_ranges = None
        self.list_of_segments_ranges_id = None
        self.lengths_of_segments = None
        self.list_of_segment_seq = None
        self.list_of_segment_db = None
        self.list_of_angles = None
        self.planar_angle = None
