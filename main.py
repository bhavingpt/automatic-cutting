import cortex
import nibabel as nib
import numpy

import struct

import os, sys, glob
import subprocess
import shlex
import utils, copy

import multiprocessing

from dotenv import load_dotenv

load_dotenv()

# regenerate the reference directory for a subject
def generate(subject, hemisphere, points):
    my_id = subject + "-" + hemisphere
    seams, walls, pts = utils.read_manual(points - 1, subject, hemisphere)

    child_dirs = next(os.walk('.'))[1]
    valid_subdirs = []

    for y in child_dirs:
        if y.endswith("-" + hemisphere) and y != my_id and len(glob.glob(y + "/*.asc")) != 0:
            valid_subdirs.append(y)

    if len(valid_subdirs) == 0: # if there are no other matching hemi dirs to line up with
        utils.generate_asc_files(subject, hemisphere, seams, walls, points - 1, pts)
        return

    reference = valid_subdirs[0]
    reference_bases = []
    reference_walls = []

    # generate reference walls
    for idx in range(5):
        nums = []
        j = 0
        while True:
            if os.path.exists(reference + "/wall" + str(idx + 1) + "_" + str(j) + ".asc"):
                nums.append(j)
                j += 1
            else:
                break

        middle = nums[int((len(nums) - 1)/2)]
        with open(reference + "/wall" + str(idx + 1) + "_" + str(middle) + ".asc") as f:
           lines = [x[:-1] for x in f.readlines()]
           for line in lines:
               content = line.split(" ")
               if len(content) > 3 and content[-1] != "0.00000": # find the point that is the base
                   reference_walls.append(int(content[0])) # put that base in 'ref bases'
    
    # generate reference bases
    for idx in range(5):
        with open(reference + "/cut" + str(idx + 1) + "_0.asc") as f:
           lines = [x[:-1] for x in f.readlines()]
           for line in lines:
               content = line.split(" ")
               if len(content) > 3 and content[-1] != "0.00000": # find the point that is the base
                   reference_bases.append(int(content[0])) # put that base in 'ref bases'
    
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

        data = [0 for _ in range(len(pts))]
        data[base] = 1
        
        if os.path.exists('temp.asc'):
            os.remove('temp.asc')

        with open("temp.asc", "w+") as f:
            inds = range(len(data))
            for ind, (x, y, z), d in zip(inds, pts, data):
                f.write("%3.3d %2.5f %2.5f %2.5f %2.5f\n" % (ind, x, y, z, d))

        FNULL = open(os.devnull, 'w')
        subj_surf_dir = os.environ['SUBJECTS_DIR'] + "/" + subject + "/surf/"

        subprocess.call(["mris_convert", "-c",
               "./temp.asc",
               subj_surf_dir + hemisphere + ".white",
               "temp_converted"], stdout = FNULL, stderr = subprocess.STDOUT)

        subprocess.call(["mri_surf2surf",
               "--srcsubject", subject,
               "--srcsurfval", subj_surf_dir + hemisphere + ".temp_converted",
               "--trgsubject", reference,
               "--trgsurfval", hemisphere + ".temp_transformed",
               "--hemi", hemisphere,
               "--trg_type", "curv"], stdout = FNULL, stderr = subprocess.STDOUT)

        locations = nib.freesurfer.read_morph_data(ref_surf_dir + hemisphere + ".temp_transformed")

        base_transformed = numpy.argmax(locations)
        print('transformed to ' + str(base_transformed))

        dists = r_surf.approx_geodesic_distance(base_transformed, m=10)
        base_dists = [dists[k] for k in reference_bases]

        for i in range(len(base_dists)): # fixes nan problem
            if numpy.isnan(base_dists[i]):
                base_dists[i] = 1000
        print(base_dists)

        seam_correspondence[idx] = base_dists
        os.system("rm temp.asc")

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

        data = [0 for _ in range(len(pts))]
        data[base] = 1
        
        if os.path.exists('temp.asc'):
            os.remove('temp.asc')

        with open("temp.asc", "w+") as f:
            inds = range(len(data))
            for ind, (x, y, z), d in zip(inds, pts, data):
                f.write("%3.3d %2.5f %2.5f %2.5f %2.5f\n" % (ind, x, y, z, d))

        FNULL = open(os.devnull, 'w')
        subj_surf_dir = os.environ['SUBJECTS_DIR'] + "/" + subject + "/surf/"

        subprocess.call(["mris_convert", "-c",
               "./temp.asc",
               subj_surf_dir + hemisphere + ".white",
               "temp_converted"], stdout = FNULL, stderr = subprocess.STDOUT)

        subprocess.call(["mri_surf2surf",
               "--srcsubject", subject,
               "--srcsurfval", subj_surf_dir + hemisphere + ".temp_converted",
               "--trgsubject", reference,
               "--trgsurfval", hemisphere + ".temp_transformed",
               "--hemi", hemisphere,
               "--trg_type", "curv"], stdout = FNULL, stderr = subprocess.STDOUT)

        locations = nib.freesurfer.read_morph_data(ref_surf_dir + hemisphere + ".temp_transformed")

        base_transformed = numpy.argmax(locations)
        print('transformed to ' + str(base_transformed))

        dists = r_surf.approx_geodesic_distance(base_transformed, m=10)
        base_dists = [dists[k] for k in reference_walls]

        for i in range(len(base_dists)): # fixes nan problem
            if numpy.isnan(base_dists[i]):
                base_dists[i] = 1000
        print(base_dists)

        wall_correspondence[idx] = base_dists
        os.system("rm temp.asc")

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

    return new_seams, new_walls

            # EB01 lh is working fine (0 confusion)
            # EB01 rh is working fine
            # EB03 lh is producing fail             TODO only four walls wtf
            # EB03 rh is working fine (0 confusion)
            # EB04 lh is working fine (0 confusion)
            # EB04 rh is working fine (0 confusion)
    
    # now that seams and walls are ordered properly - we can proceed
    utils.generate_asc_files(subject, hemisphere, new_seams, new_walls, points - 1, pts)

############################################################

def calc_points(subject):
    for x in os.walk("./" + subject):
        files = [x for x in x[2] if x.endswith(".asc")]
        break

    min_val = 100
    for i in range(1, 6):
        min_val = min(len([x for x in files if x.startswith("cut" + str(i) + "_")]), min_val)
    for i in range(1, 6):
        min_val = min(len([x for x in files if x.startswith("wall" + str(i) + "_")]), min_val)

    return -1 if min_val < 3 else min_val

def parse_reference(hemi):
    for x in os.walk("."):
        subdirs = [y for y in x[1] if y.endswith("-" + hemi)]
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

def find_match(target_subject, surface, subjects, pts, target_file):
    target_surf_dir = os.environ['SUBJECTS_DIR'] + "/" + target_subject + "/surf/"
    estimates = []
    uuidc = target_file[:-4] + "_converted"
    uuidt = target_file[:-4] + "_transformed"

    FNULL = open(os.devnull, 'w')

    for subj_id in subjects:
        subj, hemisphere = subj_id.split("-")
        subj_surf_dir = os.environ['SUBJECTS_DIR'] + "/" + subj + "/surf/"
        
        with open(subj_id + "/" + target_file) as f:
            lines = [x[:-1] for x in f.readlines()]
            for line in lines:
                content = line.split(" ")
                if len(content) > 3 and content[-1] != "0.00000":
                    pass
                    #print("Input point: " + content[0])

        subprocess.call(["mris_convert", "-c",
               "./" + subj_id + "/" + target_file, 
               subj_surf_dir + hemisphere + ".white",
               uuidc], stdout = FNULL, stderr = subprocess.STDOUT)

        subprocess.call(["mri_surf2surf", 
               "--srcsubject", subj,
               "--srcsurfval", subj_surf_dir + hemisphere + "." + uuidc,
               "--trgsubject", target_subject,
               "--trgsurfval", hemisphere + "." + uuidt,
               "--hemi", hemisphere,
               "--trg_type", "curv"], stdout = FNULL, stderr = subprocess.STDOUT)

        locations = nib.freesurfer.read_morph_data(target_surf_dir + hemisphere + "." + uuidt)
        estimates.append(numpy.argmax(locations))

    distances = {}
    for estimate in estimates:
        distances[estimate] = surface.approx_geodesic_distance(estimate, m=10)

    answer = min(pts, key = lambda pt: sum([ distances[est][pt] for est in estimates ]) )
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

    inpath = patch_filepath
    outpath = inpath + ".flat"

    save_every = 20
    save_every_str = ' -w %d'%save_every

    cmd = "mris_flatten -O fiducial{save_every_str} {inpath} {outpath}".format(inpath=inpath, outpath=outpath, save_every_str=save_every_str)

    print(cmd)
    
    p = subprocess.call(shlex.split(cmd))

    return verts

############################################################

def autocut(subject, hemisphere):
    subjects, points = parse_reference(hemisphere)
    v = cortex.Vertex.empty(subject)
    hemi = v.left if hemisphere == "lh" else v.right
    pts, polys = cortex.db.get_surf(subject, "inflated", hemisphere)
    surface = cortex.polyutils.Surface(pts, polys)
    
    todos = ["cut1_", "cut2_", "cut3_", "cut4_", "cut5_",
             "wall1_", "wall2_", "wall3_", "wall4_", "wall5_"]

    transforms = []

    # calculate and add cuts and walls
    for idx, base in enumerate(todos):
        for i in range(0, points):
            transforms.append((subject, surface, subjects, pts, base + str(i) + ".asc"))

    with multiprocessing.Pool(processes=len(transforms)) as pool:
        results = pool.starmap(find_match, transforms)

    segments = []
    idx = 0
    for base in todos:
        segment = []
        for i in range(0, points):
            segment.append(results[idx])
            idx += 1
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
        for i in range(points - 1):
            path = cortex.polyutils.Surface.geodesic_path(surface, s[i], s[i+1])
            hemi[path] = num + 1

    generate_patch(surface, subject, hemisphere, pts, points, mwall_edge, seam, smore)

    cortex.webshow(v, open_browser=False)

def main():
    os.environ["PYTHONWARNINGS"] = "ignore"
    if len(sys.argv) == 3:
        autocut(sys.argv[1], sys.argv[2])
    elif len(sys.argv) == 4:
        generate(sys.argv[1], sys.argv[2], int(sys.argv[3]))

if __name__ == "__main__":
    main()
