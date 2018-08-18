# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMPlates.jl/blob/master/LICENSE

using Base.Test
using FEMPlates

X = Dict(1 => [0.0, 0.0, 0.0],
         2 => [1.0, 0.0, 0.0],
         3 => [1.0, 1.0, 0.0],
         4 => [0.0, 1.0, 0.0],
         5 => [0.0, 0.0, 0.1],
         6 => [1.0, 0.0, 0.1],
         7 => [1.0, 1.0, 0.1],
         8 => [0.0, 1.0, 0.1])

element = Element(Hex8, [1, 2, 3, 4, 5, 6, 7, 8])
update!(element, "geometry", X)
ip = (0.0, 0.0, 0.0)
time = 0.0
X1m, X2m, X3m, X4m, X1p, X2p, X3p, X4p = element("geometry", time)

p0 = 1/2*[X1p+X1m; X2p+X2m; X3p+X3m; X4p+X4m]
pmu = 1/2*[X1p-X1m; X2p-X2m; X3p-X3m; X4p-X4m]
xi, eta, mu = ip
N1 = 0.25*(1.0-xi)*(1.0-eta)
N2 = 0.25*(1.0+xi)*(1.0-eta)
N3 = 0.25*(1.0+xi)*(1.0+eta)
N4 = 0.25*(1.0-xi)*(1.0+eta)
Nm = 0.5*(1.0-mu)
Np = 0.5*(1.0+mu)
I = eye(3)
N = [I*N1 I*N2 I*N3 I*N4]

X0 = N*p0# + mu*

a0 = 1/4*(I*X1 + I*X2 + I*X3 + I*X4)
a1 = 1/4*(-I*X1 + I*X2 + I*X3 - I*X4)
a2 = 1/4*(I*X1 - I*X2 + I*X3 - I*X4)
a3 = 1/4*(-I*X1 - I*X2 + I*X3 + I*X4)
X0 = a0 + a1*xi + a2*xi*eta + a3*eta
