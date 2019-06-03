# Autocut

## Purpose

This repository was created to automatically cut and flatten brains using previously cut brains as a reference.

## Usage

To use this work for academic purposes, please reach out to Dr Alex Huth at huth@cs.utexas.edu. This work is not licensed for anything commercial.

## Setup

These steps are described below, but to get this working you should:
  - create a .env file specifying SUBJECTS_DIR and FREESURFER_HOME
  - have installed pycortex from source and change to the `fsaverage_transform` branch
  - overwrite pycortex/cortex/blender/blendlib.py with save.txt

### Reference Directories

The code works by storing past brains in reference directories at the root level - for example, the beginnings of three reference directories are included in this repo. It requires the naming scheme SJ-XH, where SJ is the subject and XH is right or left hemisphere. Each directory should also contain the output_seam.txt and output_wall.txt which specifies where Amanda's manual cuts were.
Any reference directories at the top level will be read, so for debugging if you don't want them to be visible move them to a subdir like /inactive/.

### Obtaining output_seam.txt

The file save.txt is a commented out version of `pycortex/cortex/blender/blendlib.py`, which will write out the output_seam.txt and output_wall.txt used in an already flattened brain. To generate these .txt files you should:
  - overwrite blendlib.py with save.txt
  - reinstall pycortex using `python -m setup.py install` or something similar
  - run a Blender visualization using `cortex.segment.cut_surface("SJ", "rh", "flatten_cut5_2")`
  
## Structure

This repo exposes two functions in `main.py`:
  1. `generate(subject, hemi, n)` -> populates an existing reference dir with cut files given n intermediate points along each cut
  2. `autocut(subject, hemi)`     -> goes through the autocutting process for a given subject / hemi, using available reference dirs
  
`main.py` calls `utils.py`, which contains some simple graph search algorithms. The input to a reference directory is a list of all the points that comprise the seams and the wall, respectively. These are separated into five segments by first separating the seams into disconnected groups, and then partitioning the walls.

## Generating reference dirs

The function `utils.read_manual()` outputs the seams and walls, but these could be in any order. It's important that reference dirs have cuts in the same order, since the way autocutting is done is by averaging estimates of where intermediate points are in each reference brain.

The code checks to see if there is already a reference dir and if so transforms points from one brain to the other to align the cuts. This was initially done by using `mri_surf2surf`, but now this is done using `cortex.freesurfer.get_mri_surf2surf_matrix`, which can be found in an updated branch of pycortex. This uses a transform matrix, which is generated once and written to the reference directory for future use.

## Autocutting

The autocutting code works by parsing all the reference directories and figuring out which to use (must have the same number of intermediate points). It then writes out a patch file using `generate_patch()`, and flattens it using `mri_flatten`, which is a two-hour process.

If the flattening fails, there are some helper methods in `vis.py` to look at the patch file and debug. If the flattening succeeds (!), the code in `distort.py` will generate graphs of areal and metric distortion.

## Visualizing result

The output flattened file is `$SUBJECTS_DIR/LW/surf/rh.autocut.patch.flat`, but it should be copied to `$SUBJECTS_DIR/LW/surf/rh.autocut.flat.patch.3d` on whichever machine you want to visualize it on. After this:
  - run `cortex.freesurfer.import_flat("LW", "autocut", hemis=['rh'], sname="LW2")`
  - if you get path errors, just nuke `overlays.xvg` in the pycortex db and try it again
  - Then run either 
    - `cortex.webshow(cortex.Vertex.empty("LW2"), recache=True)` or 
    - `cortex.webshow(cortex.Vertex.empty("LW2"), recache=True)`
