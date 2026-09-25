import numpy as np #Mathematical calculations
from ase.io import read #Read atomic data
from ase.data import covalent_radii, atomic_numbers #Collect atomic data
from io_utils import get_structure
import os

#Calculate maximum covalent radius based on element
def max_radius(atoms):
    symbols = atoms.get_chemical_symbols()
    radii = [covalent_radii[atomic_numbers[s]] for s in symbols]
    rmax = 2*max(radii)
    return rmax

#Create a vector that indicates the direction of the shift in 
def build_shift(atoms,get_int):

    from plotting import view_cleanup #Create a temporary figure(s)
    
    shift = np.zeros(3) #Initiate variable with zeroes
    d = max_radius(atoms) #Use the maximum radius/diameter

    #Input a message to subsequently preview the structure
    print("\nPreviewing structure.")
    view_cleanup(atoms)

    #Print a message to present the maximum radius, especially in case of changes to the relative atom
    print(f"\nManual shift using max_radius = {d:.3f} Å")

    #Input shifts to move the system based on the number of Ni atoms. Sign is important here as well. Negative sign means left/down shift. Positive sign means right/up shift.
    nx = get_int("Shift in x (multiples of diameter, N.B. negative = left, positive = right): ")
    ny = get_int("Shift in y (multiples of diameter, N.B. negative = down, positive = up): ")

    #Shifts multiplied here
    shift[0] = nx*d
    shift[1] = ny*d

    return shift #Return shift

def apply_shift(atoms, shift):

    atoms_translate = atoms.copy() #Cannot return atoms without changing variable name, otherwise the original structure will be retained
    atoms_translate.translate(shift) #Shift here
    atoms_translate.wrap() #Wrap back into unit cell

    #Return translated atoms
    return atoms_translate

#Find and load a reference structure. CONTCAR is preferred over POSCAR, but either valid file is usable.
def get_reference_structure(folder, repeat=(1,1,1)):

    try:
        path = get_structure(folder, os.listdir(folder))
        atoms = read(path)
        if repeat != (1,1,1):
            atoms = atoms.repeat(repeat)
        return atoms

    except Exception:
        print(f"Could not load reference structure from: {folder}")
        return None


def align_to_reference(atoms, reference_atoms, reference_symbols):

    # Make a copy so the original structure is not modified.
    aligned = atoms.copy()
    # Get fractional coordinates for periodic distance matching.
    current = aligned.get_scaled_positions()
    reference = reference_atoms.get_scaled_positions()

    for symbol in reference_symbols:

        # Find atoms of the selected reference element in both structures.
        current_indices = [i for i, s in enumerate(aligned.symbols) if s == symbol]
        reference_indices = [i for i, s in enumerate(reference_atoms.symbols) if s == symbol]

        # Keep track of reference atoms that have not yet been matched, just in case.
        unused = reference_indices.copy()

        for current_i in current_indices:

            # Find the closest reference atom, accounting for periodic boundaries.
            best_i = min(unused,key=lambda i: np.linalg.norm((current[current_i] - reference[i])[:2]- np.round((current[current_i] - reference[i])[:2])))
            # Replace the current atom's position with the matched reference position.
            aligned.positions[current_i] = reference_atoms.positions[best_i]
            unused.remove(best_i)

    return aligned

def align_to_final(atoms, final_atoms, reference_symbols):

    # Make a copy so the original structure is not modified.
    aligned = atoms.copy()

    # Find the selected reference atoms in both structures.
    current_indices = [
        i for i, s in enumerate(aligned.symbols)
        if s in reference_symbols
    ]
    final_indices = [
        i for i, s in enumerate(final_atoms.symbols)
        if s in reference_symbols
    ]
    
    # Replace the current reference-atom positions with the final-state positions.
    for current_i, final_i in zip(current_indices, final_indices):
        aligned.positions[current_i] = final_atoms.positions[final_i]

    return aligned