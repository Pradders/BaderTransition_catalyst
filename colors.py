#Use this file to access functions that alter colors

#Import functions for coloring
from ase.data.colors import jmol_colors
from ase.data import atomic_numbers, chemical_symbols
from matplotlib.colors import TwoSlopeNorm, is_color_like #Used to create color gradient

import numpy as np #Mathematical calculations

#Collect atom colors
def get_atom_colors(atoms,element_colors=None):

    # Validate user-defined colours, if provided
    if element_colors is not None:

        if not isinstance(element_colors, dict):
            raise ValueError("element_colors must be a dictionary.")

        # Check each supplied element and its corresponding colour
        for symbol, colour in element_colors.items():

            # Check that the supplied element symbol is valid
            if symbol not in chemical_symbols[1:]:
                raise ValueError(f"Invalid chemical element symbol in element_colors: '{symbol}'.")

            # Check that the supplied colour is valid
            if not is_color_like(colour):
                raise ValueError(f"Invalid colour for element '{symbol}': {colour}.")
            
    #Initialise colour array/dictionary
    colors = {}
    #Go through each image to collect elements
    for i, atom in enumerate(atoms):
        # user-defined colors (highest priority)
        if element_colors is not None and atom.symbol in element_colors:
            colors[i] = element_colors[atom.symbol]
        # fallback to Jmol defaults
        else:
            colors[i] = jmol_colors[atomic_numbers[atom.symbol]]

    #print(colors) #Print colour list if desired
    
    return colors #Output colors

#Collect Bader shading colors
def get_delta_colors(delta, max_abs, cmap, tol=0.005, repeat = (1,1,1)):

    #Must multiply Bader charge by number of periodic cell repetitions
    nrep = repeat[0] * repeat[1] * repeat[2]
    delta = np.array(delta)
    delta = np.tile(delta, nrep)

    # apply tolerance, set to 0 if below threshold
    delta_plot = delta.copy()
    delta_plot[np.abs(delta_plot) < tol] = 0.0

    #Use maximum absolute value as termini for colour shading
    if len(delta_plot) == 0: #Check first that the data array is not empty
        raise ValueError("Empty delta array — nothing to plot")
    #else:
        #max_abs = np.max(np.abs(delta_plot))

    if max_abs == 0: #Check as well that the maximum values is 0, otherwise assign a miniscule number
        max_abs = 1e-12

    #Set colourmap range
    norm = TwoSlopeNorm(
        vmin=-max_abs, #Where possible, use the maximum absolute number obtained across all images as the range
        vcenter=0,
        vmax=max_abs)

    #Create colormap array
    colors = cmap(norm(delta_plot))[:, :3]

    #Return colormap and range
    return colors, norm