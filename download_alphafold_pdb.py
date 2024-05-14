import requests
import os

def download_alphafold_pdb(uniprot_id, output_dir="pdbs"):
    """Downloads the AlphaFold PDB file for a given UniProt ID.

    Args:
        uniprot_id (str): The UniProt ID of the protein.
        output_dir (str, optional): Directory to save the PDB file. Defaults to "pdbs".
    """

    base_url = "https://alphafold.ebi.ac.uk/files/AF-"

    pdb_url = f"{base_url}{uniprot_id}-F1-model_v4.pdb"  # Assuming full-length model

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    response = requests.get(pdb_url)
    if response.status_code == 200:
        pdb_file = os.path.join(output_dir, f"{uniprot_id}.pdb")
        with open(pdb_file, "wb") as f:
            f.write(response.content)
        print(f"Downloaded PDB file for {uniprot_id} to {pdb_file}")
    else:
        print(f"PDB file not found for {uniprot_id}")

    return pdb_file
