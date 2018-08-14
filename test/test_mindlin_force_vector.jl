# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMPlates.jl/blob/master/LICENSE

using Base.Test
using FEMPlates

# Create 1x1 element, assemble force vector.

X = Dict(1 => [0.0, 0.0], 2 => [1.0, 0.0],
         3 => [1.0, 1.0], 4 => [0.0, 1.0])
element = Element(Quad4, [1, 2, 3, 4])
update!(element, "geometry", X)
update!(element, "youngs modulus", 4880.0*6.0*3.0*4.0)
update!(element, "thickness", 1.0)
update!(element, "poissons ratio", 1/3)
update!(element, "distributed load", 10.0)
problem = Problem(MindlinPlate, "test plate", 3)
add_elements!(problem, [element])
time = 0.0
assemble!(problem, time)
f = full(problem.assembly.f)
@test isapprox(f[1:3:end], 2.5*[1.0, 1.0, 1.0, 1.0])
