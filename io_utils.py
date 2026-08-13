#For data acquisition

import numpy as np #Mathematical calculations
import json #JS file
import os #Operating system

#Collect CONTCAR folder location first, failing that POSCAR folder location, or failing that raise an error
def get_structure(folder, files):
    if "CONTCAR" in files:
        return os.path.join(folder, "CONTCAR")
    elif "POSCAR" in files:
        return os.path.join(folder, "POSCAR")
    else:
        raise FileNotFoundError(f"No CONTCAR or POSCAR in {folder}")

#Collect data from ACF.dat as an array
def read_bader(acf_file):

    if not os.path.isfile(acf_file): #What if the ACF.dat file is missing?
        raise FileNotFoundError(f"Missing ACF.dat: {acf_file}")

    with open(acf_file) as f: #Read file line-by-line
        lines = f.readlines()

    # find header with "CHARGE"
    try:
        header = next(l for l in lines if "CHARGE" in l.upper()) #Find header containing charge
    except StopIteration:
        raise ValueError(f"No CHARGE header found in {acf_file}")

    #Split header and index
    cols = header.split()

    #State required column headers
    required_badercols = ["CHARGE", "X", "Y", "Z"]
    #List anything missing based on the "required columns"
    missing = [col for col in required_badercols if col not in cols]
    #Check indeed if missing has any elements in it. If so, stop.
    if missing:
        raise ValueError(
        f"ACF.dat missing columns {missing} in {acf_file}"
        )

    #Try to read columns
    try:
        charge_col = cols.index("CHARGE") #Charge column
        x_col = cols.index("X") #Coordinate columns
        y_col = cols.index("Y")
        z_col = cols.index("Z")
    except ValueError as e:
        raise ValueError(f"ACF.dat format issue in {acf_file}: {e}")

    #Initialise data arrays
    charges = []
    coordinates = []

    #Sift through each line
    for line in lines:
        parts = line.split() #split lines by whitespace
        if parts and parts[0].isdigit(): #Exclude extraneous data, especially if not in columns
            coordinates.append([float(parts[x_col]),float(parts[y_col]),float(parts[z_col])]) #Append coordinates to new array
            charges.append(float(parts[charge_col])) #Append charges to new array

    return np.array(coordinates), np.array(charges) #Return new arrays

#Find initial and final folders for comparisons
def find_transition(base=os.getcwd(),initial=("ini",),final=("fin",)):
    
    structure_files = [] #Initialise array

    for root, dirs, files in os.walk(base): #Walk from base directory, i.e., where the python file is

        #Lower case for consistency
        initial = tuple(i.lower() for i in initial)
        final = tuple(f.lower() for f in final)

        #Find the initial and final states
        ini_name = next((d for d in dirs if d.lower() in initial), None)
        fin_name = next((d for d in dirs if d.lower() in final), None)

        # look for folders containing ini + fin
        if ini_name and fin_name:

            #After finding the initial and final folders, join them alongside their rooted folder directories
            ini_dir = os.path.join(root, ini_name)
            fin_dir = os.path.join(root, fin_name)

            #Collect file names in these directories
            ini_files = os.listdir(ini_dir)
            fin_files = os.listdir(fin_dir)

            #Collect POSCAR/CONTCAR directories
            ini_struct = get_structure(ini_dir, ini_files)
            fin_struct = get_structure(fin_dir, fin_files)

            #Collect ACF.dat directory
            ini_acf = os.path.join(ini_dir, "ACF.dat")
            fin_acf = os.path.join(fin_dir, "ACF.dat")

            #Root folder, relative to code location
            transition = os.path.relpath(root, base)

            #Add directories and root folder into a file for easier access
            if all([ini_struct and fin_struct and os.path.isfile(ini_acf) and os.path.isfile(fin_acf)]): #As long that all files are accessible, else skip an image(s)
                structure_files.append({
                    "transition": transition, 
                    "ini_structure": ini_struct, #Initial structure
                    "fin_structure": fin_struct, #Final structure
                    "ini_acf": ini_acf, #Initial charges
                    "fin_acf": fin_acf #Final charges
                })

            

    #Check item lists just in case to ensure that each file was indeed read
    #for item in structure_files:
    #    print(item["transition"])
    #    print("  ini:", item["ini_structure"])
    #    print("  fin:", item["fin_structure"])

    return structure_files #Output the file locations

#If Mode 3, save any constants in a json file to reuse it
def save_json(entry, filename):
    data = {
        "value": entry.tolist() if hasattr(entry, "tolist") else entry
    }
    with open(filename, "w") as f:
        json.dump(data, f)
#If Mode 3, load any constants from the saved json file to reuse it
def load_json(filename):
    try:
        with open(filename, "r") as f:
            return np.array(json.load(f)["value"])
    except FileNotFoundError:
        return None
#If Mode 3, delete the saved json file after each file is looped through
def delete_file(filename):
    try:
        os.remove(filename) #Delete if possible
    except FileNotFoundError:
        pass