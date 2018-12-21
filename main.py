import cortex
import nibabel as nib
import numpy

import os, sys
import subprocess
import utils

# regenerate the reference directory for a subject
def generate(subject, hemisphere, points):
    seams, walls, pts = utils.approximate(points - 1, subject, hemisphere)

    for x in os.walk("."):
        subdirs = [y for y in x[1] if y.endswith("-" + hemi)]
        if len(subdirs) == 0:
            utils.generate_asc_files(subject, hemisphere, seams, walls, points - 1)
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
        # TODO figure out which of reference bases it matches closest

        # TODO put that in correspondence

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

    return subjects, min(usable_points)

def find_match(target_subject, surface, subjects, points, target_file):
    target_surf_dir = os.environ['SUBJECTS_DIR'] + "/" + target_subject + "/surf/"
    estimates = []

    FNULL = open(os.devnull, 'w')

    for subj_id in subjects:
        subj, hemisphere = subj_id.split("-")
        subj_surf_dir = os.environ['SUBJECTS_DIR'] + "/" + subj + "/surf/"
        
        with open(subj_id + "/" + target_file) as f:
            lines = [x[:-1] for x in f.readlines()]
            for line in lines:
                content = line.split(" ")
                if len(content) > 3 and content[-1] != "0.00000":
                    print("Input point: " + content[0])

        subprocess.call(["mris_convert", "-c",
               "./" + subj_id + "/" + target_file, 
               subj_surf_dir + hemisphere + ".white",
               "cut_converted"], stdout = FNULL, stderr = subprocess.STDOUT)

        subprocess.call(["mri_surf2surf", 
               "--srcsubject", subj,
               "--srcsurfval", subj_surf_dir + hemisphere + ".cut_converted",
               "--trgsubject", target_subject,
               "--trgsurfval", hemisphere + ".cut_transformed",
               "--hemi", hemisphere,
               "--trg_type", "curv"], stdout = FNULL, stderr = subprocess.STDOUT)

        locations = nib.freesurfer.read_morph_data(target_surf_dir + hemisphere + ".cut_transformed")
        estimates.append(numpy.argmax(locations))

    print("    Output point: " + str(estimates[0]))
    return estimates[0]

############################################################

def autocut(subject, hemisphere):
    subjects, points = parse_reference(hemisphere)
    v = cortex.Vertex.empty(subject)
    hemi = v.left if hemisphere == "lh" else v.right
    pts, polys = cortex.db.get_surf(subject, "inflated", hemisphere)
    surface = cortex.polyutils.Surface(pts, polys)
    
    todos = ["cut1_", "cut2_", "cut3_", "cut4_", "cut5_",
             "wall1_", "wall2_", "wall3_", "wall4_", "wall5_"]

    # calculate and add cuts and walls
    for idx, base in enumerate(todos):
        print("Calculating " + base[:-1])
        segments = []
        for i in range(0, points):
            segments.append(find_match(subject, surface, subjects, points, base + str(i) + ".asc"))
        for i in range(points - 1):
            path = cortex.polyutils.Surface.geodesic_path(surface, segments[i], segments[i+1])
            hemi[path] = idx + 1

    cortex.webshow(v)

def main():
    autocut(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()
