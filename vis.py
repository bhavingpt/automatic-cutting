from main import generate
import sys
import cortex

subj = sys.argv[1]
hemi = sys.argv[2]

s, w = generate(subj, hemi, 5)

v = cortex.Vertex.empty(subj)
obj = v.left if hemi == "lh" else v.right

for i in range(5):
    for j in s[i]:
        obj[j] = i + 1
    for j in w[i]:
        obj[j] = i + 1

cortex.webshow(v, open_browser=False)
