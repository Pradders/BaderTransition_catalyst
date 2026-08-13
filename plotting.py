#Plotting functions
#When defining each of the functions, where possible or necessary, default values are set, unless a new input is allocated in the main file

#Use to collect, visualise atomic structures
from ase.visualize.plot import plot_atoms
#Use to write temporary files
from ase.io import write

#Plot atoms and use colormap where necessary
import matplotlib.pyplot as plt #Plot adsorption surfaces

import os #Operating system

#Import external functions
from colors import get_atom_colors, get_delta_colors
from layouts import create_axes, iter_axes

#Create a temporary file to check images
def view_cleanup(atoms, filename="temp_view.png", pause=True):

    #Write structure to temporary image
    write(filename, atoms,format="png",rotation="0x,0y,0z",show_unit_cell=1)

    #Display image
    img = plt.imread(filename)
    plt.figure()
    plt.imshow(img)
    plt.axis("off")
    plt.show(block=False)

    #Pause for inspection
    if pause:
        input("Press Enter to close the figure...")

    #Close figure
    plt.close()

    #Delete temp file
    try:
        os.remove(filename)
    except OSError:
        pass

#Plot all configurations, including with Bader charge difference shading gradients
def plot_bader_result(res, delta_max=1, tol=0.005, cmap=None, repeat = (1,1,1), save_dir="Bader_plots",views=None,element_colors=None,layout="mixed_left"):

    # skip empty delta
    if res["delta"] is None or len(res["delta"]) == 0:
        return
    
    #Set views to default if not given
    if views is None:
        views = [('0x,0y,0z')]

    #Select a default colormap if necessary
    if cmap is None:
        cmap = plt.cm.RdBu_r

    #Read both initial and final configurations
    atoms_ini = res["ini_structure"].copy()
    atoms_fin = res["fin_structure"].copy()
    atoms_fin_catalyst = res["fin_catalyst"].copy() #Only one is needed, since the catalyst will be the same on both sides and the transition is being modelled here.

    #For visualisation
    if repeat != (1, 1, 1):
        atoms_ini = atoms_ini.repeat(repeat)
        atoms_fin = atoms_fin.repeat(repeat)
        atoms_fin_catalyst = atoms_fin_catalyst.repeat(repeat)

    #Different rotations of images
    n_rot = len(views)

    fig, axes, mode = create_axes(layout, n_rot) #Figure arrangement

    # left plot (element colors), as set in main file
    atom_colors = get_atom_colors(atoms_fin,element_colors)
    # right plot (delta colors), use colourmap as set in main file
    slope_colors, norm = get_delta_colors(res["delta"], delta_max, cmap, tol, repeat)

    #Either horizontal or vertical
    if mode == "grid":

        for i,rotation in enumerate(views):

            if layout == "horizontal": #Horizonal array of figures

                #Initial
                plot_atoms(atoms_ini, axes[i,0], rotation=rotation,
                           show_unit_cell=0, colors=atom_colors)
                axes[i,0].set_title("Initial")

                #Final
                plot_atoms(atoms_fin, axes[i,1], rotation=rotation,
                           show_unit_cell=0, colors=atom_colors)
                axes[i,1].set_title("Final")

                #Charge transition
                plot_atoms(atoms_fin_catalyst, axes[i,2], rotation=rotation,
                           show_unit_cell=0, colors=slope_colors)
                axes[i,2].set_title("Charge transition")

            elif layout == "vertical": #Vertical array of figures

                #Initial
                plot_atoms(atoms_ini, axes[0,i], rotation=rotation,
                           show_unit_cell=0, colors=atom_colors)
                axes[0,i].set_title("Initial")

                #Final
                plot_atoms(atoms_fin, axes[1,i], rotation=rotation,
                           show_unit_cell=0, colors=atom_colors)
                axes[1,i].set_title("Final")

                #Charge transition
                plot_atoms(atoms_fin_catalyst, axes[2,i], rotation=rotation,
                           show_unit_cell=0, colors=slope_colors)
                axes[2,i].set_title("Charge transition")
    
    elif mode == "mixed": #Mixed format (horizontal, vertical)

        for i, rotation in enumerate(views):

            #Initial
            plot_atoms(atoms_ini, axes["ini"][i], rotation=rotation,
                       show_unit_cell=0, colors=atom_colors)
            axes["ini"][i].set_title("Initial")

            #Final
            plot_atoms(atoms_fin, axes["fin"][i], rotation=rotation,
                       show_unit_cell=0, colors=atom_colors)
            axes["fin"][i].set_title("Final")

        # ΔBader only drawn once (shared)
        plot_atoms(atoms_fin_catalyst, axes["delta"], rotation=views[0],
                   show_unit_cell=0, colors=slope_colors)
        axes["delta"].set_title("Charge transition")

    #Remove borders and tick marks
    for ax in iter_axes(axes):
        #ax.set_xticks([]) #Tick marks, uncheck if border should remain
        #ax.set_yticks([])
        ax.set_axis_off()

    # Set colorbar
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])

    # Generate color bar, connected to third plot
    if mode == "mixed": #Either charge plot at end 
        plt.colorbar(sm, ax=axes["delta"], label="Δ Bader charge (e)")
    else: #Or isolated charge plot
        plt.colorbar(sm, ax=axes[:,2] if layout=="horizontal" else axes[2,:],
                      label="Δ Bader charge (e)")

    # Make save directory
    os.makedirs(save_dir, exist_ok=True)

    #Use directory name to name image and remove undesirable separators
    parts = res["transition"].split(os.sep)

    #Build folder
    base_folder = os.path.join(save_dir,parts[0])
    #Use *parts if separating further
    os.makedirs(base_folder, exist_ok=True)

    #Filename
    filename = "_".join(parts) + ".png"
    #Save file to desired directory
    save_path = os.path.join(base_folder, filename)

    #Save and close figures
    plt.savefig(save_path, dpi=300, bbox_inches="tight")
    plt.close(fig)

    #Display figures if desired
    #plt.show()