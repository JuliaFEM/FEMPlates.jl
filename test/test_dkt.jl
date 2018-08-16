# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMPlates.jl/blob/master/LICENSE

using Base.Test
using FEMPlates

# Create 1x1 element, assemble stiffness matrix.

X = Dict(1 => [0.0, 0.0],
         2 => [1.0, 0.0],
         3 => [0.0, 1.0])
element = Element(Tri3, [1, 2, 3])
update!(element, "geometry", X)
update!(element, "youngs modulus", 4.0*4880.0)
update!(element, "thickness", 1.0)
update!(element, "poissons ratio", 1/3)
update!(element, "distributed load", 6.0)
problem = Problem(DKT, "test plate", 3)
add_elements!(problem, [element])
time = 0.0
assemble!(problem, time)
K = full(problem.assembly.K)
f = full(problem.assembly.f)

K_expected = [
 4880     0     0  -2440   1220  -1220  -2440  1220     0
    0   915  -305    610    305    610   -610  -610   305
    0  -305   915    610    305   -610   -610   610   305
-2440   610   610   2440      0    610      0  -610   610
 1220   305   305      0    610   -305  -1220   305   305
-1220   610  -610    610   -305    915    610  -915     0
-2440  -610  -610      0  -1220    610   2440  -610  -610
 1220  -610   610   -610    305   -915   -610   915     0
    0   305   305    610    305      0   -610     0   305]
@test isapprox(K, K_expected)

f_expected = zeros(9)
f_expected[1:3:end] = 1.0
@test isapprox(f, f_expected)
