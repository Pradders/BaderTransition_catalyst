from ase.io import read #Use to collect, visualise atomic structures
import numpy as np #Mathematical calculations
from inputs import confirm
from collections import Counter

def check_catalyst_consistency(item):

    #Atoms in each structure
    ini_atoms = read(item["ini_structure"])
    fin_atoms = read(item["fin_structure"])

    for name, atoms in [("Initial", ini_atoms), ("Final", fin_atoms)]:
        print(f"\n{name} structure")
        for i, atom in enumerate(atoms):
            x, y, z = atom.position
            print(f"{i+1:4d}  {atom.symbol:2s}  {x:10.4f} {y:10.4f} {z:10.4f}")

    # Allow the user to confirm the correspondence in atomic and elemental positions
    if not confirm(
        "\nCan you identify corresponding catalyst atoms in the two structures? (y/n): "):
        raise ValueError("Catalyst correspondence not confirmed. Analysis cancelled.")

    #Now that the systems are visually confirmed, now input the atomic indices
    catalyst_indices = sorted(set(
        i-1
        for part in input("\nEnter catalyst atom indices (e.g. 1-35,40-45,50,52): ").split(",")
        for i in (range(int(part.split("-")[0]), int(part.split("-")[1]) + 1)
                  if "-" in part else [int(part)])))

    # Check that the indices are not beyond the scope 
    if any(i < 0 or i >= len(ini_atoms) or i >= len(fin_atoms)
           for i in catalyst_indices):
        raise ValueError(
            "Catalyst atom index is outside one of the structures."
        )

    for i in catalyst_indices: # Now check that the elements are consistent between initial and final systems
        if ini_atoms[i].symbol != fin_atoms[i].symbol:
            raise ValueError(
                f"Element mismatch at atom {i+1}: "
                f"{ini_atoms[i].symbol} != {fin_atoms[i].symbol}")

    # Check positions and print the maximum displacement
    distances = [np.linalg.norm(ini_atoms[i].position - fin_atoms[i].position) for i in catalyst_indices]
    print(f"  Maximum displacement: {max(distances, default=0):.3f} Å")

    # Confirm catalyst correspondence
    if not confirm(
        "\nDo these positions correspond sufficiently for the comparison? (y/n): "):
        raise ValueError("Catalyst correspondence not confirmed. Analysis cancelled.")

    return catalyst_indices #Return the indices now that all checks have passed

#Check that Bader coordinates and POSCAR/CONTCAR coordinates are equivalent
def check_bader_alignment(atoms, acf_coords, tol=1e-3):

    #N.B. This tolerance is to not be confused with that of Bader charge difference.
    #This difference is to check that atoms in POSCAR/CONTCAR and ACF.dat are the same.
    #Here, difference should not be too large, and hence it is defaulted.

    pos = atoms.get_positions()  # Cartesian from ASE

    #Calculate difference based on absolute distance
    diff = np.linalg.norm(pos - acf_coords, axis=1)

    #Should be within threshold
    if np.max(diff) > tol:
        raise ValueError( #Else raise an error
            f"ACF.dat coordinates do not match POSCAR/CONTCAR "
            f"(max diff = {np.max(diff):.4f} Å)"
        )

#Check that the reference structure contains the requested elements and that the number of reference atoms matches the NEB structure.
def check_reference_structure(atoms, reference_atoms, reference_symbols):

    #Check that reference symbols were provided.
    if not reference_symbols:
        raise ValueError("Reference_symbols must contain at least one element symbol.")

    #Check that every requested reference element exists in the reference structure.
    reference_counts = Counter(reference_atoms.get_chemical_symbols())

    for symbol in reference_symbols:

        if reference_counts[symbol] == 0: #Should contain at least one symbol
            raise ValueError(f"Reference structure does not contain the requested element '{symbol}'.")

    #Check that the NEB structure contains the same number of reference atoms as the supplied reference structure.
    atoms_counts = Counter(atoms.get_chemical_symbols())

    for symbol in reference_symbols: #Provide data about the missing elements
        if atoms_counts[symbol] != reference_counts[symbol]:
            raise ValueError(
                f"Reference mismatch for '{symbol}': "
                f"reference contains {reference_counts[symbol]} atoms, "
                f"but the NEB structure contains {atoms_counts[symbol]}.")