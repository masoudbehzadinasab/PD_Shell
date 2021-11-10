#!/usr/bin/env python3

import mesh_functions as ms
import numpy as np

# Specify geometry and meshing info
O=np.array([-5.0,-26.0,0]) # starting point of the plate
Plate_Dim=np.array([10.0,52.0,0]) # dimensions of the plate

num_points=np.array([10,52,1]) # number of discretized points along the plate

# Call geometry generation with discretization and geometry parameters
G=ms.generate_unif_particles(O, O+Plate_Dim, num_points)

ms.save_geometry(G)
