import cortex
import nibabel as nib
import numpy
from scipy import sparse

import struct

import os, sys, glob
import subprocess
import shlex
import utils, copy

import multiprocessing

from dotenv import load_dotenv

# You need to populate SUBJECTS_DIR and FREESURFER_HOME
# these can be set in a .env file
load_dotenv()

# this method generates the reference directory for subject/hemi with n intermediate points
def generate(subject, hemisphere, points):
    my_id = subject + "-" + hemisphere
    seams, walls, pts = utils.read_manual(points - 1, subject, hemisphere)

    child_dirs = next(os.walk('.'))[1]
    valid_subdirs = []

    for y in child_dirs:
        if y.endswith("-" + hemisphere) and y != my_id and len(glob.glob(y + "/*.npy")) != 0:
            valid_subdirs.append(y)

    print(valid_subdirs)

    if len(valid_subdirs) == 0: # if there are no other matching hemi dirs to line up with
        print(seams[0])
        utils.generate_npy_files(subject, hemisphere, seams, walls, points - 1, pts)
        return

    reference = valid_subdirs[0]
    print('using ' + reference)
    reference_bases = []
    reference_walls = []

    if not os.path.exists(my_id + '/convert_' + reference.split("-")[0] + '.npz'):
        # we need to generate the transformation matrix
        matrix = cortex.freesurfer.get_mri_surf2surf_matrix(subject, hemisphere, "inflated", reference.split("-")[0])
        sparse.save_npz(my_id + '/convert_' + reference.split('-')[0] + '.npz', matrix)

    matrix = sparse.load_npz(my_id + '/convert_' + reference.split("-")[0] + '.npz')

    # generate reference walls
    for idx in range(5):
        nums = []
        j = 0
        while True:
            if os.path.exists(reference + "/wall" + str(idx + 1) + "_" + str(j) + ".npy"):
                nums.append(j)
                j += 1
            else:
                break

        middle = nums[int((len(nums) - 1)/2)]
        wall = int(numpy.load(reference + "/wall" + str(idx + 1) + "_" + str(middle) + ".npy"))
        reference_walls.append(wall) # put that base in 'ref walls'
    
    # generate reference bases
    for idx in range(5):
        base = int(numpy.load(reference + "/cut" + str(idx + 1) + "_0.npy"))
        reference_bases.append(base)

    r_pts, r_polys = cortex.db.get_surf(reference.split("-")[0], "inflated", reference.split("-")[1])
    r_surf = cortex.polyutils.Surface(r_pts, r_polys)

    print("\n\nusing " + reference + " as a reference dir")
    print("ref bases: " + str(reference_bases))

    reference = reference.split("-")[0]
    ref_surf_dir = os.environ['SUBJECTS_DIR'] + "/" + reference + "/surf/"

    seam_correspondence = dict() # will store the mapping from old base to new one
    wall_correspondence = dict()

    ##################################################################################

    for idx in range(5):
        base = seams[idx][0]
        print('\nconverting ' + str(base))

        all_pts = [0 for i in range(len(pts))]
        all_pts[base] = 1

        base_transformed = numpy.argmax(matrix * all_pts)

        print('transformed to ' + str(base_transformed))

        dists = r_surf.approx_geodesic_distance(base_transformed, m=10)
        base_dists = [dists[k] for k in reference_bases]

        for i in range(len(base_dists)): # fixes nan problem
            if numpy.isnan(base_dists[i]):
                base_dists[i] = 1000
        print(base_dists)

        seam_correspondence[idx] = base_dists

    seam_matching = dict()
    unassigned = []

    for destination in range(5):
        transforms = 0
        transformer = None
        for source in range(5):
            if source in seam_correspondence and numpy.argmin(seam_correspondence[source]) == destination:
                transforms += 1
                transformer = source

        # if easy, then assign to matching dict
        if transforms == 1:
            del seam_correspondence[transformer]
            seam_matching[transformer] = destination
        else:
            unassigned.append(destination)

    if len(unassigned) >= 3:
        print("Error: couldn't transform bases with this reference subject, cuts are overlapping.")
        exit(0)
    elif len(unassigned) == 2:
        source_a = list(seam_correspondence.keys())[0]
        source_b = list(seam_correspondence.keys())[1]

        one = seam_correspondence[source_a][unassigned[0]] + seam_correspondence[source_b][unassigned[1]]
        two = seam_correspondence[source_b][unassigned[0]] + seam_correspondence[source_a][unassigned[1]]

        seam_matching[source_a] = unassigned[0] if one < two else unassigned[1]
        seam_matching[source_b] = unassigned[1] if one < two else unassigned[0]

    ############################################################################

    print('\n\n\n\n\n\n')
    print(reference_walls)

    for idx in range(5):
        wall = walls[idx]
        base = wall[int((len(wall) - 1)/2)]
        print('\nconverting ' + str(base))

        all_pts = [0 for i in range(len(pts))]
        all_pts[base] = 1

        base_transformed = numpy.argmax(matrix * all_pts)

        print('transformed to ' + str(base_transformed))

        dists = r_surf.approx_geodesic_distance(base_transformed, m=10)
        base_dists = [dists[k] for k in reference_walls]

        for i in range(len(base_dists)): # fixes nan problem
            if numpy.isnan(base_dists[i]):
                base_dists[i] = 1000
        print(base_dists)

        wall_correspondence[idx] = base_dists

    wall_matching = dict()
    unassigned = []
    print(wall_correspondence)

    for destination in range(5):
        transforms = 0
        transformer = None
        for source in range(5):
            if source in wall_correspondence and numpy.argmin(wall_correspondence[source]) == destination:
                transforms += 1
                transformer = source

        # if easy, then assign to matching dict
        if transforms == 1:
            del wall_correspondence[transformer]
            wall_matching[transformer] = destination
        else:
            unassigned.append(destination)

    if len(unassigned) >= 3:
        print("Error: couldn't transform bases with this reference subject, walls are overlapping.")
        exit(0)
    elif len(unassigned) == 2:
        source_a = list(wall_correspondence.keys())[0]
        source_b = list(wall_correspondence.keys())[1]

        one = wall_correspondence[source_a][unassigned[0]] + wall_correspondence[source_b][unassigned[1]]
        two = wall_correspondence[source_b][unassigned[0]] + wall_correspondence[source_a][unassigned[1]]

        wall_matching[source_a] = unassigned[0] if one < two else unassigned[1]
        wall_matching[source_b] = unassigned[1] if one < two else unassigned[0]

    ############################################################################

    # reorder the seams and walls according to the dictionary
    new_seams = [[], [], [], [], []]
    new_walls = [[], [], [], [], []]

    for idx in range(5):
        new_seams[seam_matching[idx]] = seams[idx]
        new_walls[wall_matching[idx]] = walls[idx]

    for s in new_seams:
        print(s)
    print('\n')

    for w in new_walls:
        print(w)

    # now that seams and walls are ordered properly - we can proceed
    utils.generate_npy_files(subject, hemisphere, new_seams, new_walls, points - 1, pts)

############################################################

# this method finds out how many intermediate points were used in an existing ref directory
def calc_points(subject):
    for x in os.walk("./" + subject):
        files = [x for x in x[2] if x.endswith(".npy")]
        break

    min_val = 100
    for i in range(1, 6):
        min_val = min(len([x for x in files if x.startswith("cut" + str(i) + "_")]), min_val)
    for i in range(1, 6):
        min_val = min(len([x for x in files if x.startswith("wall" + str(i) + "_")]), min_val)

    return -1 if min_val < 3 else min_val

# this method reads an existing reference directory
def parse_reference(subject, hemi):
    my_id = subject + "-" + hemi
    for x in os.walk("."):
        subdirs = [y for y in x[1] if y.endswith("-" + hemi) and y != my_id]
        break

    calc = [calc_points(subject) for subject in subdirs]
    subjects = []
    usable = []

    for i in range(len(subdirs)):
        if calc[i] != -1:
            subjects.append(subdirs[i]) 
            usable.append(calc[i])

    if len(subjects) == 0:
        raise Exception("No valid references found!")

    points = max(usable, key=usable.count)
    searchdirs = [subjects[i] for i in range(len(subjects)) if usable[i] == points]

    return subjects, points

# this method converts points from one subject to another
def find_match(target_subject, surface, subjects, pts, target_file):
    estimates = []

    for subj_id in subjects:
        source_subj, source_hemi = subj_id.split("-")
        source_pts, _ = cortex.db.get_surf(source_subj, "inflated", source_hemi)

        matrix = sparse.load_npz(subj_id + '/convert_' + target_subject + '.npz')
        source_point = int(numpy.load(subj_id + '/' + target_file))

        pts = [0 for i in range(len(source_pts))]
        pts[source_point] = 1

        result = matrix * pts
        estimates.append(numpy.argmax(result))

    distances = {}
    for estimate in estimates:
        distances[estimate] = surface.approx_geodesic_distance(estimate, m=10)
        print(distances[estimate][:10])
        print('\n')

    all_points = range(len(pts))
    answer = min(all_points, key = lambda pt: sum([ distances[est][pt] for est in estimates ]) )

    print(str(answer) + " " + str(estimates))
    return answer

def generate_patch(surface, subject, hemisphere, subj_pts, intermeds, mwall_edge, seam, smore):
    all_points = set(range(len(subj_pts)))
    region_a = set()

    queue = [next(iter(all_points))]
    while len(queue) != 0:
        current = queue.pop(0)
        if current not in mwall_edge and current not in region_a:
            region_a.add(current)
            queue += list(surface.graph.neighbors(current))

    region_b = all_points - region_a

    mwall_region = min(region_a, region_b, key = len)

    # By this point, all of the vertices and edges are in proper shape

    fverts = set(range(len(subj_pts))) - mwall_region
    
    edges = mwall_edge | (smore - seam) # all points in the edge

    verts = fverts - seam
    pts = [(v, list(subj_pts[v])) for v in verts]

    # write out the patch file

    patch_filepath = os.environ['SUBJECTS_DIR'] + "/" + subject + "/surf/" + hemisphere + ".autocut.patch"
    if os.path.exists(patch_filepath):
        os.remove(patch_filepath)
    with open(patch_filepath, "wb") as f:
        f.write(struct.pack('>2i', -1, len(pts)))
        for i, pt in pts:
            if i in edges:
                f.write(struct.pack('>i3f', -i-1, *pt))
            else:
                f.write(struct.pack('>i3f', i+1, *pt))

    # The patch file for the brain has been written out - if the flattening fails,
    # you can stop here and visualize it using vis.py

    inpath = patch_filepath
    outpath = inpath + ".flat"

    save_every = 20
    save_every_str = ' -w %d'%save_every

    cmd = "mris_flatten -O fiducial{save_every_str} {inpath} {outpath}".format(inpath=inpath, outpath=outpath, save_every_str=save_every_str)

    print(cmd)
    
    p = subprocess.call(shlex.split(cmd))

    return verts

############################################################

# this method does autocutting
def autocut(subject, hemisphere):
    subjects, points = parse_reference(subject, hemisphere)
    print("Found subjects - " + str(subjects))

    for source_subject in subjects:
        if not os.path.exists(source_subject + '/convert_' + subject + '.npz'):
            # we need to generate the transformation matrix
            s_subj, s_hemi = source_subject.split("-")
            matrix = cortex.freesurfer.get_mri_surf2surf_matrix(s_subj, s_hemi, "inflated", subject)
            sparse.save_npz(source_subject + '/convert_' + subject + '.npz', matrix)

    v = cortex.Vertex.empty(subject)
    hemi = v.left if hemisphere == "lh" else v.right
    pts, polys = cortex.db.get_surf(subject, "inflated", hemisphere)
    surface = cortex.polyutils.Surface(pts, polys)
    
    todos = ["cut1_", "cut2_", "cut3_", "cut4_", "cut5_",
             "wall1_", "wall2_", "wall3_", "wall4_", "wall5_"]

    transforms = []

    segments = []
    for idx, base in enumerate(todos):
        segment = []
        for i in range(0, points):
           val = find_match(subject, surface, subjects, pts, base + str(i) + ".npy")
           segment.append(val)
        segments.append(segment)

    for i in range(6, 10): # edit one: make sure that the medial wall is continous
        segments[i][0] = segments[i - 1][-1]

    ######################
    
    mwall_edge = set()
    seam = set()
    for wall in segments[5:]:
        for i in range(points - 1):
            path = cortex.polyutils.Surface.geodesic_path(surface, wall[i], wall[i+1])
            mwall_edge.update(path)
       
    for i in range(0, 5): # edit two: make sure that seam base is on closest point on mwall
        dists = cortex.polyutils.Surface.approx_geodesic_distance(surface, segments[i][0])
        segments[i][0] = min(mwall_edge, key = lambda x: dists[x])

    for cut in segments[:5]:
        for i in range(points - 1):
            path = cortex.polyutils.Surface.geodesic_path(surface, cut[i], cut[i+1])
            seam.update(path)

    smore = set()
    for cut_point in seam:
        smore.update(surface.graph.neighbors(cut_point))
        smore.add(cut_point)
 
    ######################
    
    for num, s in enumerate(segments):
        for i in mwall_edge:
            hemi[i] = 1
        for i in smore:
            hemi[i] = 2
        for i in seam:
            hemi[i] = 3

    generate_patch(surface, subject, hemisphere, pts, points, mwall_edge, seam, smore)

def main():
    os.environ["PYTHONWARNINGS"] = "ignore"
    if len(sys.argv) == 3:
        autocut(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        generate(sys.argv[1], sys.argv[2], int(sys.argv[3]))

if __name__ == "__main__":
    main()
