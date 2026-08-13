# BaderCharge_CatalystSurface

Leverages ASE, numpy and pyplot to identify corresponding catalyst atoms between two VASP structures, calculate their differences in Bader charges, and colour the catalyst atoms according to their relative charge changes.

The code is intended for systems where the total number of atoms may differ between initial and final structures, but the catalyst atoms being compared remain the same. For example:

Ni → Ni-furfural
Ni-furfural → Ni-furfural-H

Version 1.0 (v1.0) current release.

# General procedure

## Overview

This analysis identifies corresponding catalyst atoms in the initial and final structures and calculates the Bader charge difference for these atoms only. Adsorbate and other non-catalyst atoms are excluded from the charge-transition calculation.

The catalyst atoms are manually identified after the initial and final structures are displayed. The code checks that:

1. The selected indices exist in both structures.
2. The elements at the selected indices are consistent.
3. The user confirms that the selected atoms correspond sufficiently.

The maximum displacement of the selected catalyst atoms is displayed, but identical positions are not required because structural relaxation may occur.

## Folder and file access

The atomic positions (POSCAR/CONTCAR) and Bader charges (ACF.dat) are collected from VASP calculations.

These are placed within the initial and final folders of each transition, relative to `main.py`. The default folder names are `"ini"` and `"fin"`, with additional names such as `"initial"`, `"final"`, `"is"` and `"fs"` available.

reactionORtransition/

-> ini/

->-> CONTCARorPOSCAR

->-> ACF.dat

-> fin/

->-> CONTCARorPOSCAR

->-> ACF.dat

CONTCAR/POSCAR and ACF.dat are mandatory.

N.B. The total number of atoms does **not** need to be identical between the two structures. Only the selected catalyst atoms need to correspond.

## Catalyst identification

The atomic index, element and Cartesian coordinates of every atom are displayed for both structures.

The user confirms that corresponding catalyst atoms can be identified and then enters their indices. Individual indices, ranges, or combinations can be used, e.g.:

1-35,40-45,50,52

Displayed indices begin at 1 for convenience, while the code converts them to Python's zero-based indices.

The element at each selected index is checked between the two structures. The maximum catalyst displacement is then displayed, allowing the user to decide whether the correspondence is suitable.

## Bader charge difference

The catalyst-only Bader charge difference is calculated as:

Δq = q_final - q_initial

For example:

Ni-furfural − Ni

calculates the change in Bader charge of the Ni catalyst caused by furfural adsorption.

Likewise:

Ni-furfural-H − Ni-furfural

calculates the change caused by hydrogen addition.

Only the selected catalyst atoms are included in the charge-transition plot.

The maximum absolute catalyst Bader charge difference across all analysed systems is stored and used for the colour scale. Differences below the specified tolerance (default: 0.005 e) are treated as zero.

## Projection

Atoms may extend over periodic cell boundaries, particularly adsorbates. The existing shifting functionality is therefore retained.

Four modes are available:

1. Same shift per (ini, fin) pair
2. Manual shift for EACH structure
3. Same shift for ALL structures
4. NO shift to ANY structure

Shifting is applied to the complete structures before the catalyst-only structures are extracted for plotting.

## Colour coding

The output contains:

1. Initial catalyst
2. Final catalyst
3. Catalyst Bader charge transition

The initial and final catalyst structures use elemental colours. The larger third image uses a colourmap to represent catalyst Bader charge changes.

Only the selected catalyst atoms are plotted.

## Image saving

Figures are saved to `Bader_plots`.

If the transition folders are nested within additional folders, their directory structure is retained within `Bader_plots`.

# Inputs

## Bader maximum

If an existing `bader_max.json` file is found:

Existing bader_max.json found. Load saved maximum Bader charge difference? (y/n):

The user can choose whether to reuse the stored value.

## Catalyst identification

The structures are displayed first:

Initial structure

Final structure

The user is then asked:

Can you identify corresponding catalyst atoms in the two structures? (y/n):

After confirmation:

Enter catalyst atom indices (e.g. 1-35,40-45,50,52):

The maximum displacement is displayed, followed by:

Do these positions correspond sufficiently for the comparison? (y/n):

If rejected, the analysis is cancelled.

## Key files, folders and inputs

### Usage

python main.py ([`main.py`](main.py))

### Output options

layout = "horizontal", "vertical", "split"

repeat = (1,1,1)

views = [('0x,0y,0z')]

INITIAL = ("ini","initial","is")

FINAL = ("fin","final","fs")

tol = 0.005

cmp = plt.cm.RdBu_r

save_dir = "Bader_plots"

element_colors = {"Ni": "lightgray", "C": "black"}

## External function files

[`io_utils.py`](io_utils.py) — Collects file locations and JSON files.

[`check.py`](check.py) — Checks catalyst correspondence and Bader/structure consistency.

[`analysis.py`](analysis.py) — Collects Bader charges and calculates catalyst-only charge differences.

[`plotting.py`](plotting.py) — Plots the initial catalyst, final catalyst and catalyst Bader charge transition.

[`colors.py`](colors.py) — Sets elemental colours and Bader charge colourmap.

[`layouts.py`](layouts.py) — Controls image arrangement.

[`geometry.py`](geometry.py) — Handles structural shifts.

[`inputs.py`](inputs.py) — Controls user inputs and shifting modes.

# Example

A typical application is comparison of a Ni catalyst surface before and after furfural adsorption:

Ni → Ni-furfural

and subsequent reaction steps:

Ni-furfural → Ni-furfural-H

In each case, only the corresponding Ni catalyst atoms are used for the Bader charge comparison.

The final figure contains two smaller elemental images showing the catalyst before and after the transition and one larger image showing the Bader charge changes across the catalyst surface.
