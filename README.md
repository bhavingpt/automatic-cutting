# automatic-cutting

## Purpose

This repository was created to automatically cut and flatten brains using previously cut brains as a reference. You can reach out to Bhavin with any questions (on messenger, or bhavingpt@gmail.com)!

## Setup

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

## Autocutting

## Other things
