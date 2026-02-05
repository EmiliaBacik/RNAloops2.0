import glob
import os
import subprocess
from pathlib import Path
import uuid
import shutil
import urllib.request
import json


LOG_DIR = Path("logs_from_last_update")

def update_repos():
    task = subprocess.Popen(['ls', '|', 'xargs', '-I{}', 'git', '-C', '{}', 'pull'])
    task.wait()

def run_set_update():
    print('start nonRedundantRNASetDownload')
    log_file = LOG_DIR / "log_nonRedundantRNASetDownload.txt"
    with open(log_file, "w") as log:
        subprocess.run([ 'python3', 'nonRedundantRNASetDownload/main.py'], check=True, stdout=log, stderr=log)

def run_parallel_molecule_filter():
    print('start molecule_filter')
    script = Path("nonRedundantRNASetDownload/run_molecule_filter_parallel.sh")
    if not script.exists():
        raise FileNotFoundError("Brak skryptu run_molecule_filter_parallel.sh")
    log_file = LOG_DIR / "log_molecule-filter.txt"
    with open(log_file, "w") as log:
        subprocess.run(["chmod", "+x", str(script)], check=True, stdout=log, stderr=log)
        subprocess.run(["./" + str(script)], check=True, stdout=log, stderr=log)

def run_rnapdbee():
    if not os.path.exists(Path('./dotbracket_files')):
        os.makedirs(Path('dotbracket_files'))
    for files in glob.glob('PDB_files/*.cif'):
        files = Path(files)
        if not (os.path.isfile('./dotbracket_files/' + os.path.split(str(files.with_suffix('')))[1] + '-2D-dotbracket.dbn')):
            temp_path = Path('/tmp/' + str(uuid.uuid4()))
            os.makedirs(temp_path)
            task = subprocess.Popen([ '/opt/rnapdbee-standalone-old/rnapdbee', '-i', files.absolute(), '-o', temp_path, '-a', 'DSSR'], 
stdout=subprocess.PIPE)
            task.wait()
            try:
                shutil.move(temp_path / '0' / 'strands.dbn','./dotbracket_files/' + files.stem + '-2D-dotbracket.dbn')
            except:
                f= open("log.txt","a+")
                f.write(str(files.absolute()) + "\n")
                f.close()

def dbn_cleaner():
    print('start dbn_cleaner')
    for file in glob.glob('dotbracket_files/*.dbn'):
        with open(file,"r") as dot_bracket_file:
            lines = dot_bracket_file.readlines()
            for line_number in range(len(lines)):
                if lines[line_number].startswith(">strand"):
                    dot_bracket_representation = lines[line_number + 2]
                    fragments_to_remove = [pos for pos, char in enumerate(dot_bracket_representation) if char == '-']
                    lines[line_number + 1] = "".join([char for idx, char in enumerate(lines[line_number + 1]) if idx not in fragments_to_remove])
                    lines[line_number + 2] = "".join([char for idx, char in enumerate(lines[line_number + 2]) if idx not in fragments_to_remove])
        clean_file = open(file, "w")
        clean_file.writelines(lines)
        clean_file.close()


def run_euler_angle_calculator():
    dotbracket_dir = Path("dotbracket_files")
    bash_script = Path("run_single_records.sh")
    fail_list = Path("failed_single_records.txt")
    log_file = Path("single_records.log")
    dbn_files = sorted(dotbracket_dir.glob("*.dbn"))
    if not dbn_files:
        print("Brak plików .dbn w dotbracket_files/")
        return
    with open(bash_script, "w") as f:
        f.write("#!/bin/bash\n")
        f.write("echo \"Start: $(date)\" > single_records.log\n")
        f.write("echo \"\" > failed_single_records.txt\n\n")
        for dbn in dbn_files:
            name = dbn.name
            f.write(f"echo \"Processing {name}\" >> single_records.log\n")
            f.write(
                f"python3 ./nWayJunction_release/main.py single {dbn} "
                f">> single_records.log 2>&1 || echo \"{name}\" >> failed_single_records.txt\n\n"
            )
        f.write("echo \"Finished: $(date)\" >> single_records.log\n")
    subprocess.run(["chmod", "+x", str(bash_script)], check=True)
    subprocess.run(["./" + str(bash_script)], check=True)
    subprocess.run(["python3", "./nWayJunction_release/main.py", "merge"], check=True)
    print("Zakończono batch. Sprawdź:")
    print(" - single_records.log")
    print(" - failed_single_records.txt")

def run_parallel_euler_angle_calculator():
    print('start nWayJunction')
    script = Path("nWayJunction_release/run_nWayJunction_parallel.sh")
    if not script.exists():
        raise FileNotFoundError("Brak skryptu run_nWayJunction_parallel.sh")
    log_file = LOG_DIR / "log_nWayJunction.txt"
    with open(log_file, "w") as log:
        subprocess.run(["chmod", "+x", str(script)], check=True, stdout=log, stderr=log)
        subprocess.run(["./" + str(script)], check=True, stdout=log, stderr=log)
        subprocess.run(["python3", "./nWayJunction_release/main.py", "merge"], check=True, stdout=log, stderr=log)

def run_parallel_RNApdbee():
    #task = subprocess.call(['find PDB_files/ -name \'*.cif\' > list_of_files_for_RNApdbee.txt'], shell=True, stdout=subprocess.PIPE)  
    #task.wait()
    #task = subprocess.call(['cat list_of_files_for_RNApdbee.txt | parallel --verbose --jobs 4 python3 RNApdbeeParallel/main.py'], shell=True, stdout=subprocess.PIPE)
    #task.wait()
    print('start RNApdbee')
    script = Path("RNApdbeeParallel/run_RNApdbee_parallel.sh")
    if not script.exists():
        raise FileNotFoundError("Brak skryptu run_RNApdbee_parallel.sh")
    log_file = LOG_DIR / "log_RNApdbee.txt"
    with open(log_file, "w") as log:
        subprocess.run(["chmod", "+x", str(script)], check=True, stdout=log, stderr=log)
        subprocess.run(["./" + str(script)], check=True, stdout=log, stderr=log)

def run_updater():
    print('start updater')
    log_file = LOG_DIR / "log_updater.txt"
    with open(log_file, "w") as log:
        subprocess.run(["yarn", "--cwd", "rna-loops-EH/updater", "start"], check=True, stdout=log, stderr=log)
    print('end')

def change_update_flag(value):
    url = 'http://localhost:5000/api/settings'
    data = json.dumps({
        "name": "update",
        "value": value
    }).encode('utf8')
    headers = {'content-type': 'application/json'}
    req = urllib.request.Request(url=url, data=data, method='PUT', headers=headers)
    with urllib.request.urlopen(req) as f:
        pass
    print('Status', f.status, 'Update state set to', value)

if __name__ == "__main__":
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    #update_repos()
    change_update_flag(1)
    run_set_update()
    run_parallel_molecule_filter()
    #run_rnapdbee()
    run_parallel_RNApdbee()
    dbn_cleaner()
    #run_euler_angle_calculator()
    run_parallel_euler_angle_calculator()
    run_updater()
    change_update_flag(0)
