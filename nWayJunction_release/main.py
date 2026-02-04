#!/usr/bin/python
import config
import os
from Modes.MultipleCustom import multiple_file_custom_mode, multiple_files_from_pdbee, single_files_from_pdbee
from pathlib import Path
import sys
import Commons.XML_Generator

if __name__ == "__main__":

    # Setting directories

    if not os.path.exists(Path('./output')):
        os.makedirs(Path('output'))
    if not os.path.exists(Path('./output/structures')):
        os.makedirs(Path('./output/structures'))
    if not os.path.exists(Path('./PDB_files')):
        os.makedirs(Path('PDB_files'))
    if not os.path.exists(Path('./output/single_records/')):
        os.makedirs(Path('./output/single_records/'))

    # Run main application
    if config.mode.upper() == 'MULTIPLE_CUSTOM':
        multiple_file_custom_mode(config.path_to_dotbracket_files)
    elif config.mode.upper() == 'MULTIPLE':
        if (len(sys.argv)>1):
            multiple_files_from_pdbee(sys.argv[1])
        else:
            multiple_files_from_pdbee(config.path_to_dotbracket_files)
    elif config.mode.upper() == 'CONSOLE':
        if len(sys.argv) == 0:
            print('Please pass valid arguments')
        elif sys.argv[1].upper() == 'SINGLE':
            if len(sys.argv) > 2:
                single_files_from_pdbee(sys.argv[2])
            else:
                print('No files provided!')
        elif sys.argv[1].upper() == 'MERGE':
            Commons.XML_Generator.mergeXMLs()
            Commons.XML_Generator.finalXMLs()
        else:
            print('Please pass valid arguments')
    else:
        print('Select valid mode')
