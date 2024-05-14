import imager  # Assuming your rendering functions are in a module called 'imager'
from download_alphafold_pdb import download_alphafold_pdb

def main():
    uniprot_ids = ["Q5VSL9"]  # Example UniProt IDs

    for uniprot_id in uniprot_ids:
        # 1. Download PDB
        pdb_file = download_alphafold_pdb(uniprot_id)

        # 2. Find Optimal Views
        optimal_views = imager.find_optimal_views(pdb_file)
        
        # 3. Render Images
        imager.render_molecular_surface(pdb_file, uniprot_id, optimal_views)

if __name__ == "__main__":
    main()
