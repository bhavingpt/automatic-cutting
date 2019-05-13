import cortex
from cortex.polyutils import Distortion
import matplotlib.pyplot as plt
import sys

subject = sys.argv[1]
'''
areal_dist_2 = cortex.db.get_surfinfo(subject, "distortion", dist_type="areal")
metric_dist_2 = cortex.db.get_surfinfo(subject, "distortion", dist_type="metric")
'''

# First let's load the surface and compute the distortion directly using the
# Distortion class

# load fiducial (mid-cortical) surfaces
# we're ignoring the right hemisphere surface here
# the polys (triangles) are the same for the fiducial and flat surfaces
_, (rfidpts, rpolys) = cortex.db.get_surf(subject, "fiducial")

# load flattened surfaces
_, (rflatpts, rpolys) = cortex.db.get_surf(subject, "flat")

# Create the Distortion object
dist = Distortion(rflatpts, rfidpts, rpolys)

# Compute areal distortion
# this returns an array of values for each vertex, which we will put into
# a Vertex object for plotting
areal_dist = cortex.Vertex(dist.areal, subject, vmin=-2, vmax=2)
# areal distortion is in log_2 units (e.g. -1 is half the area, 1 is double)

# Next compute metric distortion
metric_dist = cortex.Vertex(dist.metric, subject, vmin=-2, vmax=2)
# metric distortion is in mm (e.g. -1 means flatmap edge is 1 mm shorter)

cortex.quickshow(areal_dist, with_rois=False, with_labels=False)
cortex.quickshow(metric_dist, with_rois=False, with_labels=False)

# these also return Vertex objects like those we created above
plt.show()
