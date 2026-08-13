#All functions are located within external files. Many of them reference each other.
#Only three are run independently and should be called within this main file (see below)
from io_utils import find_transition #Collect file locations (POSCAR, CONTCAR, ACF.dat) for subsequent call
from analysis import collect_delta_results #Collect Bader charges and calculate charge difference between 2 states
from plotting import plot_bader_result #Plot Bader charge

import os #Operating system
import matplotlib.pyplot as plt #Plot function, called here only to account for one input variable

#Define the main function to run imported functions
def main(base, INITIAL=("ini",), FINAL=("fin",), tol=0.005, cmp=None, repeat=(1,1,1), save_dir="Bader_plots",views=None,element_colors=None,layout="split"):
    structure_files = find_transition(base,INITIAL,FINAL)
    delta_results, delta_max = collect_delta_results(structure_files, skip_errors=True)
    for res in delta_results:
        try: #Output the name of the transition here to enable the user to know which system has been processed
            parts = res["transition"].split(os.sep)
            transition_name = "_".join(parts)
            print(f"\nProcessing: {transition_name}")
        except Exception as e: #In case the file name cannot be found, pass an error message
            print(f"Error processing {res['transition']}: {e}")
            continue
        plot_bader_result(res,delta_max,tol,cmp,repeat,save_dir,views,element_colors,layout)

#Entry point/switch to run function
if __name__ == "__main__":
    base = os.getcwd()   # start from script location

    INITIAL = ("ini","initial","is") #Default is to store key files of each transition in "ini" (initial state) and "fin" (final state).
    #If there are alternative names, then please include them manually
    FINAL = ("fin","final","fs")
    tol = 0.005 #Default tolerance, i.e., if the (abs) difference in Bader charge is less than 0.005, it will treated as 0
    cmp = plt.cm.RdBu_r #Default colormap also in function, can be changed to any desired, plt already imported
    repeat = (1,1,1) #In case of adsorbate atoms extending over the unit cell, this will increase the size of periodicity
    save_dir="Bader_plots" #Save folder

    views = [ #Different rotations to view atomic/surface configurations
    ('0x,0y,0z'),     # top
    #('-90x,0y,0z'),    # side
    #('-90x,-90y,0z'), # front
    ] #Choose one of ('0x,0y,0z') #top, ('90x,0y,0z') #side, ('-90x,-90y,0z') #front, or add other rotations as desired

    element_colors = { #Define desired colors for atoms in POSCAR/CONTCAR if desired
        "Ni": "lightgray",
        "C": "black",
        #"O": "red", #O is already red by default
        #"H": "white", #H is already white by default
        }
    
    layout = "split" #also "horizontal" and "vertical"
    #"split" means that the standard ini/fin images will be above each other on the LHS of the image
    #"horizontal" arranges all images horizontally, "vertical" arranges all images vertically

    main(base, INITIAL, FINAL, tol, cmp, repeat, save_dir, views, element_colors, layout) #Start main function