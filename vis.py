from main import generate
import sys
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
