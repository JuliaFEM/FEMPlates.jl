# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/PDEAssembler.jl/blob/master/LICENSE

using JuliaFEM
using JuliaFEM.Preprocess
add_elements! = JuliaFEM.add_elements!

#field_elements, boundary_elements = get_unit_square()

meshfile = Pkg.dir("FEMPlates", "examples", "dkt", "mesh.med")
mesh = aster_read_mesh(meshfile)

using FEMPlates

field_elements = create_elements(mesh, "PLATE")
update!(field_elements, "youngs modulus", 4880.0)
update!(field_elements, "poissons ratio", 1/3)
update!(field_elements, "thickness", 0.5)
update!(field_elements, "distributed load", 10.0)
plate = Problem(DKT, "test problem", 3)
add_elements!(plate, field_elements)

boundary_elements = create_elements(mesh, "BORDER")
update!(boundary_elements, "displacement 1", 0.0)
bc = Problem(Dirichlet, "fixed", 3, "displacement")
add_elements!(bc, boundary_elements)

analysis = Analysis(Linear)
add_problems!(analysis, [plate, bc])
run!(analysis)

X = plate("geometry", 0.0)
u = plate("displacement", 0.0)
N = length(X)
x = [X[i][1] for i=1:N]
y = [X[i][2] for i=1:N]
w = [u[i][1] for i=1:N]

using Plots
surface(x, y, w)
gui()
using PyPlot
elset_name = :PLATE
node_ids = sort(collect(keys(mesh.nodes)))
node_x = [mesh.nodes[id][1] for id in node_ids]
node_y = [mesh.nodes[id][2] for id in node_ids]
element_ids = mesh.element_sets[elset_name]
element_conn = [mesh.elements[id]-1 for id in element_ids]
triplot(node_x, node_y, element_conn, "-k")
axis("equal")
axis("off")
