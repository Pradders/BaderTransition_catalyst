#Plotting functions
#When defining each of the functions, where possible or necessary, default values are set, unless a new input is allocated in the main file

#Use to collect, visualise atomic structures
from ase.visualize.plot import plot_atoms
#Use to write temporary files
from ase.io import write

#Plot atoms and use colormap where necessary
import matplotlib.pyplot as plt #Plot adsorption surfaces

import numpy as np #Mathematical analysis
import os #Operating system

#Import external functions
from colors import get_atom_colors, get_delta_colors
from layouts import create_axes, iter_axes
from geometry import align_to_reference, align_to_final

# Default labels
DEFAULT_LABELS = {
    "initial": "Initial",
    "final": "Final",
    "delta": "Δq"
}

# Default styles
DEFAULT_STYLES = {
    "heading": {
        "fontname": "Times New Roman",
        "fontsize": 24,
        "fontweight": "bold"
    },
    "tick": {
        "fontname": "Times New Roman",
        "fontsize": 18
    }
}

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
def plot_bader_result(res, delta_max=1, tol=0.005, cmap=None, repeat = (1,1,1), save_dir="Bader_plots", views=None, element_colors=None,
                      layout="mixed_left", labels=None, styles=None, reference_atoms=None, reference_symbols=("Ni",)):

    from check import check_reference_structure

    # skip empty delta
    if res["delta"] is None or len(res["delta"]) == 0:
        return
    
    #Set views to default if not given
    if views is None:
        views = [('0x,0y,0z')]

    #Select a default colormap if necessary
    if cmap is None:
        cmap = plt.cm.RdBu_r

    #Check labels. Use default if invalid or not passed.
    try:
        labels["initial"]
        labels["final"]
        labels["delta"]
    except (TypeError, KeyError):
        labels = DEFAULT_LABELS

    #Check styles. Use default if invalid or not passed.
    try:
        styles["heading"]
        styles["tick"]
    except (TypeError, KeyError):
        styles = DEFAULT_STYLES

    #Read both initial and final configurations
    atoms_ini = res["ini_structure"].copy()
    atoms_fin = res["fin_structure"].copy()
    atoms_fin_catalyst = res["fin_catalyst"].copy() #Only one is needed, since the catalyst will be the same on both sides and the transition is being modelled here.
    #The final structure is better since it is more convenient to compare to the final state.

    #Activate if reference atoms are defined
    if reference_atoms is not None:
        # Use the reference surface for the initial and final visual structures.
        for atoms in (atoms_ini, atoms_fin):
            check_reference_structure(atoms, reference_atoms, reference_symbols)
        atoms_ini = prepare_visual_structure(atoms_ini,reference_atoms,reference_symbols)
        atoms_fin = prepare_visual_structure(atoms_fin,reference_atoms,reference_symbols)

        # Align the Bader catalyst to the external reference surface.
        atoms_fin_catalyst.wrap()
        reference_atoms = reference_atoms.copy()
        reference_atoms.wrap()
        atoms_fin_catalyst = align_to_reference(atoms_fin_catalyst,reference_atoms,reference_symbols)

    #Otherwise, use the structures as they are...
    else:
        atoms_fin_catalyst.wrap()
        atoms_fin.wrap()
        atoms_ini.wrap()
        #...but the initial Ni atoms should match with those in the final image for consistency
        atoms_ini = align_to_final(atoms_ini,atoms_fin,reference_symbols)

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
                axes[i,0].set_title(labels["initial"], **styles["heading"])

                #Final
                plot_atoms(atoms_fin, axes[i,1], rotation=rotation,
                           show_unit_cell=0, colors=atom_colors)
                axes[i,1].set_title(labels["final"], **styles["heading"])

                #Change in Bader charge
                plot_atoms(atoms_fin_catalyst, axes[i,2], rotation=rotation,
                           show_unit_cell=0, colors=slope_colors)
                axes[i,2].set_title(labels["delta"], **styles["heading"])

            elif layout == "vertical": #Vertical array of figures

                #Initial
                plot_atoms(atoms_ini, axes[0,i], rotation=rotation,
                           show_unit_cell=0, colors=atom_colors)
                axes[0,i].set_title(labels["initial"], **styles["heading"])

                #Final
                plot_atoms(atoms_fin, axes[1,i], rotation=rotation,
                           show_unit_cell=0, colors=atom_colors)
                axes[1,i].set_title(labels["final"], **styles["heading"])

                #Change in Bader charge
                plot_atoms(atoms_fin_catalyst, axes[2,i], rotation=rotation,
                           show_unit_cell=0, colors=slope_colors)
                axes[2,i].set_title(labels["delta"], **styles["heading"])
    
    elif mode == "mixed": #Mixed format (horizontal, vertical)

        for i, rotation in enumerate(views):

            #Initial
            plot_atoms(atoms_ini, axes["ini"][i], rotation=rotation,
                       show_unit_cell=0, colors=atom_colors)
            axes["ini"][i].set_title(labels["initial"], **styles["heading"])

            #Final
            plot_atoms(atoms_fin, axes["fin"][i], rotation=rotation,
                       show_unit_cell=0, colors=atom_colors)
            axes["fin"][i].set_title(labels["final"], **styles["heading"])

        # ΔBader only drawn once (shared)
        plot_atoms(atoms_fin_catalyst, axes["delta"], rotation=views[0],
                   show_unit_cell=0, colors=slope_colors)
        axes["delta"].set_title(labels["delta"], **styles["heading"])

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
        cbar = plt.colorbar(sm, ax=axes["delta"])
        cbar.set_label(label=f"{labels['delta']} (e)",**styles["heading"])
    else: #Or isolated charge plot
        cbar = plt.colorbar(sm, ax=axes[:,2] if layout=="horizontal" else axes[2,:])
        cbar.set_label(label=f"{labels['delta']} (e)",**styles["heading"])

    for tick in cbar.ax.get_yticklabels():
        plt.setp(tick, **styles["tick"])
    
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

#Prepare a visual copy of an NEB image so that the selected reference atoms use the same periodic representation as the reference structure.
def prepare_visual_structure(atoms, reference_atoms, reference_symbols):

    visual_reference = reference_atoms.copy()

    #Get fractional coordinates inside the unit cell.
    current_scaled = atoms.get_scaled_positions(wrap=False)
    reference_scaled = reference_atoms.get_scaled_positions(wrap=True)

    cell = atoms.get_cell()

    #Keep track of the occurrence number of each element.
    #This allows the reference structure to contain more than
    #one type of substrate atom.
    reference_indices = {}
    current_indices = {}

    for symbol in reference_symbols:

        reference_indices[symbol] = [
            i for i, atom_symbol
            in enumerate(reference_atoms.get_chemical_symbols())
            if atom_symbol == symbol]

        current_indices[symbol] = [
            i for i, atom_symbol
            in enumerate(atoms.get_chemical_symbols())
            if atom_symbol == symbol]

    # Find one common periodic translation for the whole structure.
    best_shift = np.zeros(3)
    best_distance = np.inf   

    for x_shift in (-1, 0, 1):

        for y_shift in (-1, 0, 1):

            total_distance = 0.0

            for symbol in reference_symbols:

                for reference_index, current_index in zip(
                    reference_indices[symbol],
                    current_indices[symbol]):

                    # Current atom with a common periodic translation.
                    candidate = current_scaled[current_index].copy()
                    candidate[0] += x_shift
                    candidate[1] += y_shift

                    # Difference from the corresponding reference atom.
                    difference = candidate - reference_scaled[reference_index]

                    # Convert to Cartesian distance.
                    distance = np.linalg.norm(difference @ cell)

                    total_distance += distance**2

            # Keep the common translation giving the smallest
            # total distance for all reference atoms.
            if total_distance < best_distance:

                best_distance = total_distance
                best_shift = np.array([x_shift, y_shift, 0.0])


    # Apply the SAME periodic translation to the entire structure.
    reference_scaled += best_shift
    visual_reference.set_cell(cell)
    visual_reference.set_scaled_positions(reference_scaled)

    # Collect only the atoms that are NOT part of the reference substrate.
    non_reference_indices = [
        i for i, symbol in enumerate(atoms.get_chemical_symbols())
        if symbol not in reference_symbols]

    visual_adsorbates = atoms[non_reference_indices]

    # Combine reference substrate with the non-reference atoms.
    visual_atoms = visual_reference + visual_adsorbates

    return visual_atoms