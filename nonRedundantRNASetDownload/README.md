#README

##Description

The purpose of this tool is to get the current non-redundant set of RNAs or identify which structures have been added/removed since the last update.

The workflow consists of:

1. Downloading all PDB ID of current RNA structures from PBD database (https://www.rcsb.org/)
2. Creating the initial_set or files_to_update_file:
* If RNA_SETS folder contains only one file, that means that there is no initial set and it needs to be created. The script creates init_set.txt file that contains all structures PDB ID from the latest non-redundant set.
* If the RNA_SETS folder contains more than one file, it means that the user has already established an initial set and needs to identify files that need to be added/removed. The script creates files_to_update.txt file which contains structures that need to be added (thair PDB is is
preceded by '+' sign)or removed (PDB ID preceded by '-' sign).
3. Saving PDB ID that failed to download in the logs_from_last_update/failed_nonRedundantRNASetDownload.txt

##Run

To run type:

python main.py in the console
