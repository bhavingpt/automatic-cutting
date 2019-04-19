from main import generate
import sys, os
import cortex

def hi(subj, hemi):
    s, w = generate(subj, hemi, 5)

    v = cortex.Vertex.empty(subj)
    obj = v.left if hemi == "lh" else v.right

    pts, polys = cortex.db.get_surf(subj, "inflated", hemi)
    surface = cortex.polyutils.Surface(pts, polys)

    for i in range(5):
        seam = s[i]
        wall = w[i]

        for j in range(len(seam) - 1):
            path = cortex.polyutils.Surface.geodesic_path(surface, seam[j], seam[j+1])
            obj[path] = i + 1

        for j in range(len(wall) - 1):
            path = cortex.polyutils.Surface.geodesic_path(surface, wall[j], wall[j+1])
            obj[path] = i + 1

    return v

def com(subj, hemi, nums):
    v = cortex.Vertex.empty(subj)
    obj = v.left if hemi == "lh" else v.right

    obj[nums[0]] = 1
    for n in nums[1:]:
        obj[n] = 2

    return v

if __name__ == "__main__":
    # grab patch file from adele and vis it
    command_one = ["scp", "-P", "26452", "bhavin@129.116.157.223:/usr/local/freesurfer/subjects/LW/surf/rh.autocut.patch", "."]
    command_two = ["mv", "rh.autocut.patch", "$SUBJECTS_DIR/LW/surf/rh.autocut.patch.3d"]

    os.system(" ".join(command_one))
    os.system(" ".join(command_two))

    cortex.freesurfer.show_surf("LW", "rh", "inflated", "autocut")
