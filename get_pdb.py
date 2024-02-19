##############################################################################
import requests

def download_pdb_file(pdb_id, save_path):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(url)
    if(response.status_code == 200):
        with open(save_path, 'wb') as f:
            f.write(response.content)
        print(f"PDB file '{pdb_id}.pdb' downloaded successfully.")
    else:
        print(f"Failed to download PDB file for {pdb_id}. Status code: {response.status_code}")

with open("pdb_id.dat") as f:
    lines = f.readlines()
for line in lines:
    terms = line.split()
    if len(terms) != 3:
        continue
    pdb_id = terms[0]
    if (pdb_id == 'pdb_id'):
        continue
    save_path = f"data/pdb/{pdb_id}.pdb"
    download_pdb_file(pdb_id, save_path)

##############################################################################

