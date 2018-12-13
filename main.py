import cortex
import os, sys

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

############################################################

def autocut(subject, hemisphere):
    subjects, points = parse_reference(hemisphere)
    v = cortex.Vertex.empty(subject)

    # calculate and add the cuts
    for i in range(1, 6):
        pass

    # calculate and add the walls
    for i in range(1, 6):
        pass


def main():
    autocut(sys.argv[1], sys.argv[2])

if __name__ == "__main__":
    main()
