import imager
from download_alphafold_pdb import download_alphafold_pdb
import pandas as pd
import sys
from joblib import Parallel, delayed

def process_protein(uniprot_id):
    try:
        pdb_file = download_alphafold_pdb(uniprot_id)
        optimal_views = imager.find_optimal_views(pdb_file)
        imager.render_molecular_surface(pdb_file, uniprot_id, optimal_views)
    except Exception as e:
        print(f"Error processing {uniprot_id}: {e}")  # Log errors for debugging

def main():
    proteome_file = pd.read_excel(sys.argv[1])
    uniprot_ids = proteome_file['Entry'].tolist()

    # Parallelization
    num_jobs = -1  # Use all available CPU cores (adjust if needed)
    Parallel(n_jobs=num_jobs)(delayed(process_protein)(uniprot_id) for uniprot_id in uniprot_ids)

if __name__ == "__main__":
    main()
