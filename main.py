import cortex
import nibabel as nib
import numpy

import struct

import os, sys
import subprocess
import shlex
import utils, copy

import multiprocessing

from dotenv import load_dotenv

load_dotenv()

# regenerate the reference directory for a subject
def generate(subject, hemisphere, points):
    my_id = subject + "-" + hemisphere
    seams, walls, pts = utils.approximate(points - 1, subject, hemisphere)

    for x in os.walk("."):
        subdirs = [y for y in x[1] if y.endswith("-" + hemisphere) and y != my_id]
        if len(subdirs) == 0:
            utils.generate_asc_files(subject, hemisphere, seams, walls, points - 1, pts)
            return

    reference = subdirs[0]
    reference_bases = []

    for idx in range(5):
        with open(reference + "/cut" + str(idx) + "_0.asc") as f:
           lines = [x[:-1] for x in f.readlines()]
           for line in lines:
               content = line.split(" ")
               if len(content) > 3 and content[-1] != "0.00000":
                   reference_bases.append(int(content[0]))
    
    r_pts, r_polys = cortex.db.get_surf(reference.split("-")[0], "inflated", reference.split("-")[1])
    r_surf = cortex.polyutils.Surface(r_pts, r_polys)

    reference = reference.split("-")[0]
    ref_surf_dir = os.environ['SUBJECTS_DIR'] + "/" + reference + "/surf/"

    correspondence = dict() # will store the mapping from old base to new one

    for idx in range(5):
        base = seam[idx]

        data = [0 for _ in range(len(pts))]
        data[base] = 1
        with open("temp.asc", "w+") as f:
            inds = range(len(data))
            for ind, (x, y, x), d in zip(inds, pts, data):
                f.write("%3.3d %2.5f %2.5f %2.5f %2.5f\n" % (ind, x, y, z, d))

        FNULL = open(os.devnull, 'w')
        subj_surf_dir = os.environ['SUBJECTS_DIR'] + "/" + subject + "/surf/"

        subprocess.call(["mris_convert", "-c",
               "temp.asc",
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

        dists = r_surf.approx_geodesic_distance(base_transformed)
        base_dists = [dists[k] for k in reference_bases]

        correspondence[idx] = numpy.argmin(base_dists)
        os.system("rm temp.asc")

    # reorder the seams and walls according to the dictionary
    new_seams = [[], [], [], [], []]
    new_walls = [[], [], [], [], []]

    for idx in range(5):
        new_seams[correspondence[idx]] = seams[idx]
        new_walls[correspondence[idx]] = walls[idx]
         
    utils.generate_asc_files(subject, hemisphere, new_seams, new_walls, points - 1)

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

def find_match(target_subject, surface, subjects, points, target_file):
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

    #print("    Output point: " + str(estimates[0]))
    return estimates[0]

def generate_patch(surface, subject, hemisphere, subj_pts, intermeds, cuts, walls):
    mwall_edge = set()
    seam = set()

    for cut in cuts:
        for i in range(intermeds - 1):
            path = cortex.polyutils.Surface.geodesic_path(surface, cut[i], cut[i+1])
            seam.update(path)

    for wall in walls:
        for i in range(intermeds - 1):
            path = cortex.polyutils.Surface.geodesic_path(surface, wall[i], wall[i+1])
            mwall_edge.update(path)

    smore = set()
    for cut_point in seam:
        smore.update(surface.graph.neighbors(cut_point))
        smore.add(cut_point)

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

    save_every = 1
    save_every_str = ' -w %d'%save_every

    cmd = "mris_flatten -O fiducial{save_every_str} {inpath} {outpath}".format(inpath=inpath, outpath=outpath, save_every_str=save_every_str)

    print(cmd)
    
    subprocess.check_call(shlex.split(cmd))

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
            transforms.append((subject, surface, subjects, points, base + str(i) + ".asc"))

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

    for i in range(6, 10): # add the previous wall's end to the beginning
        segments[i][0] = segments[i - 1][-1]

    for num, s in enumerate(segments):
        if num > 4:
            print(s)
        for i in range(points - 1):
            path = cortex.polyutils.Surface.geodesic_path(surface, s[i], s[i+1])
            hemi[path] = num + 1

    generate_patch(surface, subject, hemisphere, pts, points, segments[:5], segments[5:])

    cortex.webshow(v, open_browser=False)

def main():
    autocut(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()
