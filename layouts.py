#Plot atoms and use colormap where necessary
import matplotlib.pyplot as plt #Plot adsorption surfaces
import matplotlib.gridspec as gridspec #Graphing in grid format

import numpy as np #Mathematical calculations

def create_axes(layout, n_rot):
    if layout == "horizontal": #Horizontal
        fig, axes = plt.subplots(n_rot, 3, figsize=(12, 6*n_rot), squeeze=False)
        return fig, axes, "grid" #return also a keyword

    elif layout == "vertical": #Vertical
        fig, axes = plt.subplots(3, n_rot, figsize=(6*n_rot, 10), squeeze=False)
        return fig, axes, "grid" #return also a keyword

    elif layout == "split": #Split
        fig = plt.figure(figsize=(10, 6*n_rot))
        gs = gridspec.GridSpec(2*n_rot, 2, width_ratios=[1, 1])

        #Pre-define axes arrays
        axes = {
                "ini": [],
                "fin": [],
                "delta": None
            }

        #Automatically add initial and final subplots on LHS
        for i in range(n_rot):
                axes["ini"].append(fig.add_subplot(gs[2*i, 0]))
                axes["fin"].append(fig.add_subplot(gs[2*i+1, 0]))

        #Add the charge subplot alone on the RHS
        axes["delta"] = fig.add_subplot(gs[:, 1])  # spans full right side

        #return also a keyword
        return fig, axes, "mixed"

    else:
        raise ValueError(f"Unknown layout: {layout}")

def iter_axes(axes): #generator function,  
    if isinstance(axes, dict): #Check for a dictionary
        for v in axes.values():
            if isinstance(v, list): #Convert to a list
                for ax in v:
                    yield ax #faster output
            else:
                yield v
    else:
        for ax in np.ravel(axes):
            yield ax #1D list