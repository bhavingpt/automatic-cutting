import cortex
import nibabel
import numpy

import os, sys
import subprocess

# regenerate the reference directory for a subject
def generate(subject, hemisphere, points):
    pass

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
    usable_points = []

    for i in range(len(subdirs)):
        if calc[i] != -1:
            subjects.append(subdirs[i]) 
            usable_points.append(calc[i])

    if len(subjects) == 0:
        raise Exception("No valid references found!")

    return subjects, min(usable_points)

def find_match(target_subject, surface, subjects, points, target_file):
    target_surf_dir = os.environ['SUBJECTS_DIR'] + "/" + target_subject + "/surf/"
    estimates = []

    FNULL = open(os.devnull, 'w')

    for subj_id in subjects:
        subj, hemisphere = subj_id.split("-")
        subj_surf_dir = os.environ['SUBJECTS_DIR'] + "/" + subj + "/surf/"

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

        locations = nibabel.freesurfer.read_morph_data(target_surf_dir + hemisphere + ".cut_transformed")
        estimates.append(numpy.argmax(locations))

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
        segments = []
        for i in range(1, points + 1):
            segments.append(find_match(subject, surface, subjects, points, base + str(i) + ".asc"))
        for i in range(points - 1):
            path = cortex.polyutils.Surface.geodesic_path(surface, segments[i], segments[i+1])
            hemi[path] = idx + 1

    cortex.webshow(v)

def main():
    autocut(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()
