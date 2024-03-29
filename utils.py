import cortex
import random
import os, glob
import numpy as np
import networkx as nx

from collections import deque

obj = None
surface = None

def create(direc):
    global obj
    print("Creating the vertex...")
    if not os.path.isdir(direc):
        print("Files for this subject/hemi not found")
        exit(0)
    subj = direc.split("-")[0]
    v = cortex.Vertex.empty(subj)
    obj = v.left if "lh" in direc else v.right

    with open(direc + "/output_wall.txt") as f:
        pts = [int(x[:-1]) for x in f.readlines()]
    obj[pts] = 1

    #with open(direc + "/output_smore.txt") as f:
    #    pts = [int(x[:-1]) for x in f.readlines()]
    #obj[pts] = 3

    with open(direc + "/output_seam.txt") as f:
        pts = [int(x[:-1]) for x in f.readlines()]
    obj[pts] = 2

    return v

def separate_seams(graph, all_cuts):
    seams = []
    used = set()
    up_next = list()
    current_list = -1

    for candidate in all_cuts:
        if candidate not in used:
            used.add(candidate)
            current_list += 1
            seams.append([candidate])
            up_next += list(graph.neighbors(candidate))

        while len(up_next) > 0:
            x = up_next.pop()
            if x not in used and x in all_cuts:
                used.add(x)
                seams[current_list].append(x)
                up_next += list(graph.neighbors(x))

    return seams

def separate_wall(graph, medial_wall, bases):
    global surface

    dists = dict()
    for base in bases:
        dists[base] = cortex.polyutils.Surface.approx_geodesic_distance(surface, base)

    paths = []

    used = set()

    next_base = bases[0]
    first_base = bases[0]

    while True:
        if next_base == None:
            if len(bases) == 0:
                break
            next_base = bases[0]
        bases.remove(next_base)

        prev_dict = dict()
        prev = None
        up_next = [next_base]
        path = []

        quit = False
        while not quit:
            if len(up_next) == 0:
                path.append(first_base)
                next_base = None
                break

            x = up_next.pop()
            if x in used:
                continue

            prev_dict[x] = prev
            prev = x
            used.add(x)

            for y in graph.neighbors(x):
                if y not in used and y not in bases and y in medial_wall:
                    up_next.append(y)

            for base in bases:
                if dists[base][x] < 3:
                    next_base = base 
                    prev_dict[next_base] = x
                    quit = True

        while prev is not None:
            path.append(prev)
            prev = prev_dict[prev]
        paths.append(path)

    return paths

def generate_bases(seams, medial_wall):
    global surface
    dists = cortex.polyutils.Surface.approx_geodesic_distance(surface, medial_wall)
    return [min(seam, key = lambda x: dists[x]) for seam in seams]

def get_ends(direc, segments):
    global surface

    with open(direc + "/output_wall.txt") as f:
        medial_wall = [int(x[:-1]) for x in f.readlines()]

    with open(direc + "/output_seam.txt") as f:
        all_cuts = [int(x[:-1]) for x in f.readlines()]

    seams = separate_seams(surface.graph, all_cuts)
    seam_landmarks = []
    for i in range(len(seams)):
        print("    Finding seam " + str(i) + "'s intermediate points...")
        landmarks = []
        seam = seams[i]

        ordered_seam = []
        neighbors = dict()
        for point in seam:
            neighbors[point] = []
            for end in list(surface.graph.neighbors(point)):
                if end in seam:
                    neighbors[point].append(end)

        ends = {k: v for k, v, in neighbors.items() if len(v) == 1}
        current = list(ends.keys())[0]
        used = set()

        while current != None:
            ordered_seam.append(current)
            used.add(current)
            queue = set(neighbors[current]) - used
            current = min(queue, key = lambda x: len(neighbors[x])) if len(queue) > 0 else None

        pieces = np.array_split(ordered_seam, segments)
        for piece in pieces:
            landmarks.append(piece[0])
        landmarks.append(pieces[-1][-1])
        seam_landmarks.append(landmarks)

    bases = generate_bases(seam_landmarks, medial_wall)

    for descriptor in seam_landmarks:
        if descriptor[0] not in bases:
            descriptor.reverse()

    walls = separate_wall(surface.graph, medial_wall, bases)
    wall_landmarks = []
    for i in range(len(walls)):
        print("    Finding wall " + str(i) + "'s intermediate points...")

        landmarks = []
        ordered_wall = walls[i]
        if len(ordered_wall) < segments:
            continue

        pieces = [list(x) for x in np.array_split(ordered_wall, segments)]
        for piece in pieces:
            landmarks.append(piece[0])
        landmarks.append(pieces[-1][-1])
        wall_landmarks.append(landmarks)

    for descriptor in wall_landmarks:
        descriptor.reverse()

    return seam_landmarks, wall_landmarks

def read_manual(segments, subject, hemi, style="inflated"):
    global obj
    global surface

    direc = subject + '-' + hemi
    v = create(direc)
    pts, polys = cortex.db.get_surf(subject, style, hemi)
    surface = cortex.polyutils.Surface(pts, polys)

    seams, walls = get_ends(direc, segments)

    return seams, walls, pts

    # visualize the cuts
    for j in range(segments):
        path = cortex.polyutils.Surface.geodesic_path(surface, seam[j], seam[j+1])
        obj[path] = 3 + i

    # visualize the walls
    for j in range(segments):
        path = cortex.polyutils.Surface.geodesic_path(surface, wall[j], wall[j+1])
        obj[path] = 3 + i + len(seams)

    cortex.webshow(v)

def generate_npy_files(subject, hemi, seams, walls, segments, pts):
    direc = subject + '-' + hemi

    os.chdir(direc)
    os.system("rm -rf *.npy")

    # handle the cuts
    for i, seam in enumerate(seams):
        # write out the cuts asc files
        for j in range(0, segments + 1):
            filename = "cut" + str(i + 1) + "_" + str(j) + ".npy"
            np.save(filename, np.asarray(seam[j]))

    # handle the walls
    for i, wall in enumerate(walls):
        # write out the walls asc files
        for j in range(0, segments + 1):
            filename = "wall" + str(i + 1) + "_" + str(j) + ".npy"
            np.save(filename, np.asarray(wall[j]))

    os.chdir("..")
