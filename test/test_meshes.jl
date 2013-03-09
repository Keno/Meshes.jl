using Meshes

# 1. Produce a level set function that is a noisy version of the distance from
# the origin (such that level sets are noisy spheres).
#
# The noise should exercise marching tetrahedra's ability to produce a water-
# tight surface in all cases (unlike standard marching cubes).
#
N = 10
sigma = 2.0
distance = Float32[ sqrt(float32(i*i+j*j+k*k)) for i = -N:N, j = -N:N, k = -N:N ]
distance = distance + sigma*rand(2*N+1,2*N+1,2*N+1)

# 2. Extract an isosurface.
#
lambda = N-2*sigma # isovalue
tic()
msh = isosurface(distance,lambda)
toc()

# 3. Export the mesh to a ply.
#
# The mesh can be visualized, e.g., in MeshLab (http://meshlab.sourceforge.net/).
#   
exportToPly(msh,"noisy_sphere.ply")

