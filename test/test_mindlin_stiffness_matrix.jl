# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMPlates.jl/blob/master/LICENSE

using Base.Test
using FEMPlates

# Create 1x1 element, assemble stiffness matrix.

X = Dict(1 => [0.0, 0.0], 2 => [1.0, 0.0],
         3 => [1.0, 1.0], 4 => [0.0, 1.0])
element = Element(Quad4, [1, 2, 3, 4])
update!(element, "geometry", X)
update!(element, "youngs modulus", 4880.0*6.0*3.0*4.0)
update!(element, "thickness", 1.0)
update!(element, "poissons ratio", 1/3)
problem = Problem(MindlinPlate, "test plate", 3)
add_elements!(problem, [element])
time = 0.0
assemble!(problem, time)
K = full(problem.assembly.K)
#display(K)
fixed_dofs = [1, 2, 3, 10, 11, 12]
free_dofs = setdiff(1:12, fixed_dofs)
f = zeros(12)
f[[4, 7]] = 10.0
u = zeros(f)
u[free_dofs] = K[free_dofs, free_dofs] \ f[free_dofs]
#display(u)
u_expected = [
0.0
0.0
0.0
0.000306921
-0.000249543
-0.0000458344
0.000306921
-0.000249543
0.0000458344
0.0
0.0
0.0
]
@test isapprox(u, u_expected,atol=1.0e-8)
