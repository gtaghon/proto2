import sys
sys.path.append('/usr/local/Cellar/pymol/3.0.0/libexec/lib/python3.12/site-packages')
from pymol import cmd
import numpy as np

def find_optimal_views(pdb_file):
    """Calculates 4 orthogonal views with minimal surface overlap.

    Args:
        pdb_file (str): Path to the PDB file.

    Returns:
        list: A list of tuples, each representing a view as (x, y, z) in degrees.
    """

    from Bio.PDB import PDBParser

    # Load the PDB structure
    parser = PDBParser()
    structure = parser.get_structure("protein", pdb_file)

    # Extract atom coordinates
    coordinates = np.array([atom.get_coord() for atom in structure.get_atoms()])

    # Calculate principal component analysis (PCA)
    mean_coords = np.mean(coordinates, axis=0)
    centered_coords = coordinates - mean_coords
    _, _, eigenvectors = np.linalg.svd(centered_coords)

    # Create initial orthogonal views (top/bottom, front/back) AS TUPLES
    views = [
        (np.degrees(np.arctan2(eigenvectors[0][1], eigenvectors[0][0])), 0, 0),  # Convert to tuple
        (np.degrees(np.arctan2(eigenvectors[1][1], eigenvectors[1][0])), 0, 0),  # Convert to tuple
    ]

    # Find a third view that maximizes surface exposure with minimal overlap
    best_view = None
    max_dist = -np.inf  # Initialize to negative infinity
    for angle in range(0, 360, 45):  # Test at 45-degree intervals for efficiency
        new_view = (views[0][0], angle, 0)  # Use first view's x-rotation
        dist = np.linalg.norm(np.cross(eigenvectors[0], rotate_vector(eigenvectors[1], new_view)))
        if dist > max_dist:
            max_dist = dist
            best_view = new_view
    views.append(best_view)

    # Add an opposite view for the third view
    views.append((best_view[0] + 180, best_view[1], best_view[2]))

    # Format views as list of 3-tuples (x, y, z) and add the 90-degree Y rotation
    return [(view[0], 90, view[2]) for view in views]

def rotate_vector(vector, rotation):
    """Helper function to rotate a 3D vector."""
    x, y, z = rotation
    rotation_matrix = np.array([
        [np.cos(np.radians(y))*np.cos(np.radians(z)), -np.cos(np.radians(y))*np.sin(np.radians(z)), np.sin(np.radians(y))],
        [np.sin(np.radians(x))*np.sin(np.radians(y))*np.cos(np.radians(z)) + np.cos(np.radians(x))*np.sin(np.radians(z)), -np.sin(np.radians(x))*np.sin(np.radians(y))*np.sin(np.radians(z)) + np.cos(np.radians(x))*np.cos(np.radians(z)), -np.sin(np.radians(x))*np.cos(np.radians(y))],
        [-np.cos(np.radians(x))*np.sin(np.radians(y))*np.cos(np.radians(z)) + np.sin(np.radians(x))*np.sin(np.radians(z)), np.cos(np.radians(x))*np.sin(np.radians(y))*np.sin(np.radians(z)) + np.sin(np.radians(x))*np.cos(np.radians(z)), np.cos(np.radians(x))*np.cos(np.radians(y))]
    ])
    return np.dot(rotation_matrix, vector)

def yrb(selection='all'):
	""" 
	A script to highlight hydrophobicity and charge on protein surfaces
	DHS065 Hagemans et al YRB script
	created by Dominique Hagemans and Ianthe A.E.M. van Belzen, July 2015
	Rudiger group CPC, Utrecht University

	yellow: C, CH, CH2, CH3 groups that are not bound to N or O groups.
	red: negatively charged atoms
	blue: positively charged atoms
	grey: backbone, polar groups and remaining atoms (darkened from 0.95 to 0.75)
	"""
	
	cmd.remove("hydro")
	cmd.set_color('yellow',[0.950,0.78,0.0])
	cmd.set_color('grey',[0.75,0.75,0.75])
	cmd.set_color('red',[1.0,0.4,0.4])	
	cmd.set_color('blue',[0.2,0.5,0.8])	
	
	mapping = {}
	mapping['arg'] = [ ('NE,NH2,NH1', 'blue'), ('CD,CZ', 'grey'), ('CG', 'yellow') ]
	mapping['asn'] = [ ('CG,OD1,ND2', 'grey') ]
	mapping['asp'] = [ ('CG', 'grey'), ('OD2,OD1', 'red')  ]
	mapping['cys'] = [ ('SG', 'grey') ]	
	mapping['gln'] = [ ('CG', 'yellow'), ('CD,OE1,NE2', 'grey') ]
	mapping['glu'] = [ ('CG', 'yellow'), ('CD', 'grey'), ('OE1,OE2', 'red') ]
	mapping['his'] = [ ('CG,CD2,ND1,NE2,CE1', 'grey') ]	
	mapping['ile'] = [ ('CG1,CG2,CD1', 'yellow') ]
	mapping['leu'] = [ ('CG,CD1,CD2', 'yellow') ]
	mapping['lys'] = [ ('CG,CD', 'yellow'), ('CE', 'grey'), ('NZ', 'blue') ]
	mapping['met'] = [ ('CG,CE', 'yellow'), ('SD', 'grey') ]
	mapping['phe'] = [ ('CG,CD1,CE1,CZ,CE2,CD2', 'yellow') ]
	mapping['pro'] = [ ('CG', 'yellow'), ('CD', 'grey') ]
	mapping['ser'] = [ ('CB,OG', 'grey') ]
	mapping['thr'] = [ ('CB,OG1', 'grey'), ('CG2', 'yellow') ]
	mapping['trp'] = [ ('CG,CD2,CZ2,CH2,CZ3,CE3', 'yellow'), ('CD1,NE1,CE2', 'grey') ]
	mapping['tyr'] = [ ('CG,CE1,CD1,CE2,CD2', 'yellow'), ('CZ,OH', 'grey') ]
	mapping['val'] = [ ('CG1,CG2', 'yellow') ]

	obj_list = cmd.get_names('objects')
	for obj in obj_list:
		if (obj == selection or selection == 'all'):
			cmd.color('grey','(n. N,C,CA,O and ' + obj + ')')
			cmd.color('yellow','(n. CB and ' + obj + ')')
			
			for key in mapping:
				for (atom, color) in mapping[key]:
					cmd.color(color, '( n. ' + atom + ' and r. ' + key + ' and ' + obj + ' )')

def render_molecular_surface(pdb_file, output_prefix, views, surface_type="surface"):
    """Renders molecular surface images from multiple views.

    Args:
        pdb_file (str): Path to the input PDB file.
        output_prefix (str): Prefix for output image filenames.
        views (list): List of tuples, each representing a view as (x, y, z) in degrees.
        surface_type (str, optional): Type of structure to render. Defaults to "surface".
    """
    
    cmd.load(pdb_file)
    cmd.show(surface_type)  

    # Customize the appearance (optional)
    cmd.set("transparency", 0.4)
    cmd.set("orthoscopic", "off")
    cmd.set("ray_trace_mode", 0)
    yrb()
    # cmd.set("bg_rgb", [1, 1, 1]) 

    for i, view in enumerate(views):
        # Set the viewpoint (corrected 18-element vector)
        print(i)
        print(view)
        cmd.reset() # Initial default view

        # Apply rotation based on the view
        cmd.rotate("x", view[0])
        cmd.rotate("y", view[1])
        cmd.rotate("z", view[2])
		
        # Zoom to include the entire model in frame
        cmd.zoom(complete=1)

        # Render and save the image
        output_image = f"{output_prefix}_view{i+1}.png"
        cmd.draw(width=512, height=512, antialias=0)  # Adjust dimensions and resolution as needed
        cmd.png(output_image)
        print(f"Image saved to {output_image}")

# Example usage
pdb_file = "AFmodel.pdb"  # Replace with your actual PDB file path
output_prefix = "protein_surface" 

optimal_views = find_optimal_views(pdb_file)
render_molecular_surface(pdb_file, output_prefix, optimal_views)
