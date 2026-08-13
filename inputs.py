from geometry import build_shift, apply_shift #Atomic shifting functions
from plotting import view_cleanup #Create a temporary figure(s)
import numpy as np #Mathematical calculations
from io_utils import save_json, load_json

_current_mode = None #Global variable in selecting mode of imaging

#Obtain integers for certain prompts. If not an integer, then repeat.
def get_int(prompt):
    while True:
        try:
            return int(input(prompt))
        except ValueError:
            print("Invalid input. Please enter an integer (e.g., -1, 0, 1).")

#Choose whether to apply shifts to all or some images
def choose_mode():

    global _current_mode #Global variable
    if _current_mode is not None:
        return _current_mode

    #Choose a figure collection method
    while True:
        print("\nSelect shift mode:")
        print("1: SAME shift per (ini, fin) pair") #Mode 1: same shift for each (ini, fin) image pair
        print("2: MANUAL shift for EACH image") #Mode 2: manual shift for each individual image
        print("3: SAME shift for ALL images") #Mode 3: same shift for all images
        print("4: NO shift to ANY image") #Mode 4: no applied shift to any individual image

        #Which method? Select an integer among 1, 2, and 3.
        mode = get_int("Enter mode (1/2/3/4): ")
        #Ask if correct, else restart loop
        if mode in [1, 2, 3, 4]:
            confirm_mode = input(f"Confirm mode {mode}? (y/n): ").lower() #Confirmation
            if confirm_mode == "y":
                _current_mode = mode #Keep selected method constant throughout program
                return mode

        print("Invalid mode, try again.\n")

#Used to accept a result
def confirm(prompt="Accept this result/outcome? (y/n): "):
    while True:
        ans = input(prompt).lower()
        if ans in ["y", "n"]:
            return ans == "y"
        print("Please enter 'y' or 'n'.")



#Depending on the mode selected, apply appropriate imaging collection method
def process_structures(atoms_ini, atoms_fin):

    mode = choose_mode() #Which mode?

    if mode == 1: #Mode 1: same shift for each (ini, fin) image pair

        while True:
            shift = build_shift(atoms_ini,get_int) #Prepare shift
            test_ini = apply_shift(atoms_ini, shift) #Apply shift to initial image
            test_fin = apply_shift(atoms_fin, shift) #Apply shift to initial image

            view_cleanup(test_fin) #Check images temporarily if they look OK, especially after input shift constant factor

            if confirm(): #Accept result or not
                break

        return test_ini, test_fin #Return shifted results

    elif mode == 2: #Mode 2: manual shift for each individual image

        results = [] #Initialise
        for i in [atoms_ini, atoms_fin]: #Loop through each dataset
            while True:
                shift = build_shift(i,get_int) #Prepare shift
                test = apply_shift(i,shift) #Apply shift to initial/final image
                view_cleanup(test) #Check images temporarily if they look OK, especially after input shift constant factor

                if confirm(): #Accept result or not
                    results.append(test) #Collect data
                    break

        test_ini = results[0] #Same output variables as other modes
        test_fin = results[1]

        return test_ini, test_fin #Return shifted results

    elif mode == 3: #Mode 3: same shift for each individual image

        shift = load_json("shift.json") #Load constant shift

        if shift is None: #If no shift file, save one first and then reuse it

            while True:
                shift = build_shift(atoms_ini,get_int) #Prepare shift
                test_ini = apply_shift(atoms_ini,shift) #Apply shift to initial image
                test_fin = apply_shift(atoms_fin,shift) #Apply shift to final image

                view_cleanup(test_fin) #Check images temporarily if they look OK, especially after input shift constant factor
                
                if confirm(): #Accept result or not
                    save_json(shift,"shift.json") #Save shift json file
                    break

        else:
            test_ini = apply_shift(atoms_ini, shift) #Same output variables as other modes
            test_fin = apply_shift(atoms_fin, shift)

        return test_ini, test_fin #Return shifted results
    
    elif mode == 4: #Mode 4: no applied shift to any individual image

        return atoms_ini, atoms_fin #Leave input unchanged
