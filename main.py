import imager  # Assuming your rendering functions are in a module called 'imager'
from download_alphafold_pdb import download_alphafold_pdb
import pandas as pd
import sys

def main():
    proteome_file = pd.read_excel(sys.argv[1])
    uniprot_ids = proteome_file['Entry'].tolist()
    # uniprot_ids = ["H6WB16"]  # Example UniProt IDs

    for uniprot_id in uniprot_ids:
        # 1. Download PDB
        pdb_file = download_alphafold_pdb(uniprot_id)

        # 2. Find Optimal Views
        optimal_views = imager.find_optimal_views(pdb_file)
        
        # 3. Render Images
        imager.render_molecular_surface(pdb_file, uniprot_id, optimal_views)

if __name__ == "__main__":
    main()
