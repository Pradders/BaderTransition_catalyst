#Data analysis
#When defining each of the functions, where possible or necessary, default values are set, unless a new input is allocated in the main file

#Use to collect, visualise atomic structures
from ase.io import read

import numpy as np #Mathematical calculations

#Import external functions
from io_utils import read_bader, delete_file, load_json, save_json
from check import check_catalyst_consistency, check_bader_alignment
from inputs import process_structures, confirm

import os #Operating system

#Calculate Bader charge difference
def collect_delta_results(structure_files, skip_errors=False):

    #Initialise the array containing the Bader charge differences
    delta_results = []

    #Initialise a variable to determine the maximum range of Bader charge
    delta_max = 0.0

    if os.path.exists("bader_max.json"): #In case an existing bader_max.json file can be found

        if confirm("\nExisting bader_max.json found. Load saved maximum Bader charge difference? (y/n): "):
            delta_max = float(load_json("bader_max.json")) #If desired, load bader_max file

    for item in structure_files: #Read file list

        try:

            # read structures
            ini_atoms = read(item["ini_structure"])
            fin_atoms = read(item["fin_structure"])

            #Check in case files are empty
            if len(ini_atoms) == 0 or len(fin_atoms) == 0:
                raise ValueError(f"Empty POSCAR/CONTCAR structure in {item['transition']}")

            ini_coords, ini_q = read_bader(item["ini_acf"]) #Read initial bader file and extract coordinates and charges
            fin_coords, fin_q = read_bader(item["fin_acf"]) #Read final bader file and extract coordinates and charges

            #Check in case files are empty
            if len(ini_q) == 0 or len(fin_q) == 0:
                raise ValueError(f"Empty ACF.dat file in {item['transition']}")

            #Coordinate alignment checks
            check_bader_alignment(ini_atoms, ini_coords)
            check_bader_alignment(fin_atoms, fin_coords)

            #Identify and check catalyst
            catalyst_indices = check_catalyst_consistency(item)

            #Check catalyst Bader data
            if any(i >= len(ini_q) or i >= len(fin_q) for i in catalyst_indices):
                raise ValueError("Catalyst atom index is outside the Bader data.")

            #Calculate Bader charge difference
            delta_q = np.round(np.array(fin_q)[catalyst_indices] - np.array(ini_q)[catalyst_indices], 3)

            #Calculate maximum charge difference
            local_max = np.max(np.abs(delta_q))
            #Check if the local maximum is higher than the current maximum. If so, update delta_max
            if local_max > delta_max:
                delta_max = local_max

            #Unwrap both initial and final structures
            ini_atoms_update, fin_atoms_update = process_structures(ini_atoms, fin_atoms)

            #Extract catalyst structures
            ini_catalyst = ini_atoms_update[catalyst_indices]
            fin_catalyst = fin_atoms_update[catalyst_indices]

            #Collect all data and store in array. Should plot initial and final structures to observe both atomic movement and charge differences more easily.
            delta_results.append({
                "transition": item["transition"],
                "ini_structure": ini_atoms_update,
                "fin_structure": fin_atoms_update,
                "ini_catalyst": ini_catalyst,
                "fin_catalyst": fin_catalyst,
                "catalyst_indices": catalyst_indices,
                "delta": delta_q
            })

        except Exception as e:
            if skip_errors: #Only if true, set to false otherwise
                print(e)
                continue
            else:
                raise #Raise error if necessary        

    #If the shift json file is created, then delete it after use
    if os.path.exists("shift.json"):
        print("Existing shift.json detected. Delete shift.json?")
        if confirm():
            delete_file("shift.json")
            print("shift.json deleted.")

    #Print if desired
    #print(delta_results)

    #Save updated maximum for future reuse
    save_json(delta_max, "bader_max.json")

    #Output this array for future use
    return delta_results, delta_max