from impress import FineScaleMesh as mesh
from packs.stokes_brinkman_3d import stokes_brinkman as stokes

file='./mesh/9x9x9.msh'
print("Preprocessing mesh")
M = mesh(file)
print("Preprocessing Mesh Finished")
print("Applying Stokes-Brinkmann")#
stokes.stokes_solver(M)
print("Applying stokes_brinkman finished")
