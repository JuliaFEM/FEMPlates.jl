# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMPlates.jl/blob/master/LICENSE

"""
    DKT - Discrete Kirchhoff Triangle element

# References

- Batoz, Jean‐Louis, Klaus‐JÜRgen Bathe, and Lee‐Wing Ho.
  "A study of three‐node triangular plate bending elements."
  International Journal for Numerical Methods in Engineering 15.12 (1980): 1771-1812.
- Lucena Neto, Eliseu, et al.
  "An Explicit Consistent Geometric Stiffness Matrix for the DKT Element."
  Latin American Journal of Solids and Structures 14.4 (2017): 613-628.
"""
struct DKT <: FieldProblem
    geometric_stiffness :: Bool
end

function DKT()
    geometric_stiffness = true
    return DKT(geometric_stiffness)
end

function FEMBase.get_unknown_field_name(::Problem{DKT})
    return "displacement"
end

"""
        basis_factory(::Problem{DKT}, element::Element{Tri3})

Return basis functions and their partial derivatives for DKT.
"""
function basis_factory(::Problem{DKT}, element::Element{Tri3})
    X = element("geometry", time)
    x1, x2, x3 = (X[i][1] for i in 1:3)
    y1, y2, y3 = (X[i][2] for i in 1:3)

    # for ij = 23, 31, 12:
    # x_ij = x_i - x_j
    # y_ij = y_i - y_j
    # L_ij^2 = x_ij^2 - y_ij^2
    x23 = x2-x3
    x31 = x3-x1
    x12 = x1-x2
    y23 = y2-y3
    y31 = y3-y1
    y12 = y1-y2
    L23 = (x23^2 + y23^2)
    L31 = (x31^2 + y31^2)
    L12 = (x12^2 + y12^2)

    # for k = 4,5,6 and ij=23,31,12:

    # a_k = -x_ij/L_ij^2
    # b_k = 3*x_ij*y_ij / (4*L_ij^2)
    # c_k = (x_ij^2 - 2*y_ij^2) / 4*L_ij^2
    # d_k = -y_ij / L_ij^2
    # e_k = (y_ij^2 - 2*x_ij^2) / 4*L_ij^2
    a4 = -x23/L23
    a5 = -x31/L31
    a6 = -x12/L12
    b4 = 3*x23*y23/(4*L23)
    b5 = 3*x31*y31/(4*L31)
    b6 = 3*x12*y12/(4*L12)
    c4 = (x23^2 - 2*y23^2)/(4*L23)
    c5 = (x31^2 - 2*y31^2)/(4*L31)
    c6 = (x12^2 - 2*y12^2)/(4*L12)
    d4 = -y23/L23
    d5 = -y31/L31
    d6 = -y12/L12
    e4 = (y23^2 - 2*x23^2)/(4*L23)
    e5 = (y31^2 - 2*x31^2)/(4*L31)
    e6 = (y12^2 - 2*x12^2)/(4*L12)

    # p_k = -6*x_ij/l_ij^2
    # q_k = 3*x_ij*y_ij/L_ij^2
    # r_k = 3*y_ij^2/L_ij^2
    # t_k = -6*y_ij/L_ij^2
    p4 = -6*x23/L23
    p5 = -6*x31/L31
    p6 = -6*x12/L12
    q4 =  3*x23*y23/L23
    q5 =  3*x31*y31/L31
    q6 =  3*x12*y12/L12
    r4 =  3*y23^2/L23
    r5 =  3*y31^2/L31
    r6 =  3*y12^2/L12
    t4 = -6*y23/L23
    t5 = -6*y31/L31
    t6 = -6*y12/L12

    # println("(a4, a5, a6) = ($a4, $a5, $a6)")
    # println("(b4, b5, b6) = ($b4, $b5, $b6)")
    # println("(c4, c5, c6) = ($c4, $c5, $c6)")
    # println("(d4, d5, d6) = ($d4, $d5, $d6)")
    # println("(e4, e5, e6) = ($e4, $e5, $e6)")
    # println("(p4, p5, p6) = ($p4, $p5, $p6)")
    # println("(q4, q5, q6) = ($q4, $q5, $q6)")
    # println("(r4, r5, r6) = ($r4, $r5, $r6)")
    # println("(t4, t5, t6) = ($t4, $t5, $t6)")

    function N(xi)
        xi, eta = ip.coords
        N1 = 2.0*(1.0-xi[1]-xi[2])*(0.5-xi[1]-xi[2])
        N2 = xi[1]*(2.0*xi[1]-1.0)
        N3 = xi[2]*(2.0*xi[2]-1.0)
        N4 = 4.0*xi[1]*xi[2]
        N5 = 4.0*xi[2]*(1.0-xi[1]-xi[2])
        N6 = 4.0*xi[1]*(1.0-xi[1]-xi[2])
        return [N1 N2 N3 N4 N5 N6]
    end

    function Hx(xi)
        N1, N2, N3, N4, N5, N6 = N(xi)
        Hx1 = 1.5*(a6*N6 - a5*N5)
        Hx2 = b5*N5 + b6*N6
        Hx3 = N1 - c5*N5 - c6*N6
        Hx4 = 1.5*(a4*N4 - a6*N6)
        Hx5 = b4*N4 + b6*N6
        Hx6 = N2 - c4*N4 - c6*N6
        Hx7 = 1.5*(a5*N5 - a4*N4)
        Hx8 = b4*N4 + b5*N5
        Hx9 = N3 - c4*N4 - c5*N5
        return [Hx1 Hx2 Hx3 Hx4 Hx5 Hx6 Hx7 Hx8 Hx9]
    end

    function Hy(xi)
        N1, N2, N3, N4, N5, N6 = N(xi)
        Hy1 = 1.5*(d6*N6 - d5*N5)
        Hy2 = -N1 + e5*N5 + e6*N6
        Hy3 = -b5*N5 - b6*N6
        Hy4 = 1.5*(d4*N4 - d6*N6)
        Hy5 = -N2 + e4*N4 + e6*N6
        Hy6 = -b4*N4 - b6*N6
        Hy7 = 1.5*(d5*N5 - d4*N4)
        Hy8 = -N3 + e4*N4 + e5*N5
        Hy9 = -b4*N4 - b5*N5
        return [Hy1 Hy2 Hy3 Hy4 Hy5 Hy6 Hy7 Hy8 Hy9]
    end

    function dHxdxi(ip)
        xi, eta = ip
        d1 =  p6*(1.0-2.0*xi) + (p5-p6)*eta
        d2 =  q6*(1.0-2.0*xi) - (q5+q6)*eta
        d3 = -4.0 + 6.0*(xi+eta) + r6*(1.0-2.0*xi) - eta*(r5+r6)
        d4 = -p6*(1.0-2.0*xi) + eta*(p4+p6)
        d5 =  q6*(1.0-2.0*xi) - eta*(q6-q4)
        d6 = -2.0 + 6.0*xi + r6*(1.0-2.0*xi) + eta*(r4-r6)
        d7 = -eta*(p5+p4)
        d8 =  eta*(q4-q5)
        d9 = -eta*(r5-r4)
        return [d1 d2 d3 d4 d5 d6 d7 d8 d9]
    end

    function dHydxi(ip)
        xi, eta = ip
        d1 =  t6*(1.0-2.0*xi) + eta*(t5-t6)
        d2 =  1.0 + r6*(1.0-2.0*xi) - eta*(r5+r6)
        d3 = -q6*(1.0-2.0*xi) + eta*(q5+q6)
        d4 = -t6*(1.0-2.0*xi) + eta*(t4+t6)
        d5 = -1.0 + r6*(1.0-2.0*xi) + eta*(r4-r6)
        d6 = -q6*(1.0-2*xi) - eta*(q4-q6)
        d7 = -eta*(t4+t5)
        d8 =  eta*(r4-r5)
        d9 = -eta*(q4-q5)
        return [d1 d2 d3 d4 d5 d6 d7 d8 d9]
    end

    function dHxdeta(ip)
        xi, eta = ip
        d1 = -p5*(1.0-2.0*eta) - xi*(p6-p5)
        d2 =  q5*(1.0-2.0*eta) - xi*(q5+q6)
        d3 = -4.0 + 6.0*(xi+eta) + r5*(1.0-2.0*eta) - xi*(r5+r6)
        d4 =  xi*(p4+p6)
        d5 =  xi*(q4-q6)
        d6 = -xi*(r6-r4)
        d7 =  p5*(1.0-2.0*eta) - xi*(p4+p5)
        d8 =  q5*(1.0-2.0*eta) + xi*(q4-q5)
        d9 = -2.0 + 6.0*eta + r5*(1.0-2.0*eta+xi*(r4-r5))
        return [d1 d2 d3 d4 d5 d6 d7 d8 d9]
    end

    function dHydeta(ip)
        xi, eta = ip
        d1 = -t5*(1.0-2.0*eta) - xi*(t6-t5)
        d2 =  1.0 + r5*(1.0-2.0*eta) - xi*(r5+r6)
        d3 = -q5*(1.0-2.0*eta) + xi*(q5+q6)
        d4 =  xi*(t4+t6)
        d5 =  xi*(r4-r6)
        d6 = -xi*(q4-q6)
        d7 =  t5*(1.0-2.0*eta) - xi*(t4+t5)
        d8 = -1.0 + r5*(1.0-2.0*eta) + xi*(r4-r5)
        d9 = -q5*(1.0-2.0*eta) - xi*(q4-q5)
        return [d1 d2 d3 d4 d5 d6 d7 d8 d9]
    end

    return Hx, Hy, dHxdxi, dHxdeta, dHydxi, dHydeta

end

function FEMBase.assemble_elements!(problem::Problem{DKT},
                                    assembly::Assembly,
                                    elements::Vector{Element{Tri3}},
                                    time::Float64)

    ndofs = 3*3
    Ke = zeros(ndofs, ndofs)
    fe = zeros(ndofs)

    for element in elements

        fill!(Ke, 0.0)
        fill!(fe, 0.0)

        Hx, Hy, dHxdxi, dHxdeta, dHydxi, dHydeta = basis_factory(problem, element)

        for ip in get_integration_points(element)

            detJ = element(ip, time, Val{:detJ})
            E = element("youngs modulus", ip, time)
            t = element("thickness", ip, time)
            nu = element("poissons ratio", ip, time)

            D = E*t^3/(12*(1-nu^2)) * [
                1.0  nu  0.0
                nu  1.0  0.0
                0.0 0.0  0.5*(1-nu)]

            B = 1.0/(x31*y12 - x12*y31) * [
                 y31*dHxdxi(ip) + y12*dHxdeta(ip)
                -x31*dHydeta(ip) - x12*dHydeta(ip)
                -x31*dHxdxi(ip) - x12*dHxdeta(ip) + y31*dHydxi(ip) + y12*dHydeta(ip)]

            Ke += ip.weight * B'*D*B * detJ

            # if problem.properties.geometric_stiffness
            #     # Assemble Kg
            # end

            if haskey(element, "distributed load")
                N = element(ip, time)
                p = element("distributed load", ip, time)
                fe[1:3:end] += ip.weight * p * N' * detJ
            end

        end

        gdofs = get_gdofs(problem, element)
        add!(assembly.K, gdofs, gdofs, Ke)
        add!(assembly.f, gdofs, fe)

    end

    return nothing

end
