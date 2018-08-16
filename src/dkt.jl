# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMPlates.jl/blob/master/LICENSE

"""
    DKT - Discrete Kirchhoff Triangle element

# References

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

function FEMBase.assemble_elements!(problem::Problem{DKT},
                                    assembly::Assembly,
                                    elements::Vector{Element{Tri3}},
                                    time::Float64)

    for element in elements

        ndofs = 3*3
        Ke = zeros(ndofs, ndofs)
        fe = zeros(ndofs)
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

        # for k = 4,5,6 and ij=23,31,12, respectively
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

        p4 = -6*x23/L23
        p5 = -6*x31/L31
        p6 = -6*x12/L12
        t4 = -6*y23/L23
        t5 = -6*y31/L31
        t6 = -6*y12/L12
        q4 =  3*x23*y23/L23
        q5 =  3*x31*y31/L31
        q6 =  3*x12*y12/L12
        r4 =  3*y23^2/L23
        r5 =  3*y31^2/L31
        r6 =  3*y12^2/L12

        # println("(a4, a5, a6) = ($a4, $a5, $a6)")
        # println("(b4, b5, b6) = ($b4, $b5, $b6)")
        # println("(c4, c5, c6) = ($c4, $c5, $c6)")
        # println("(d4, d5, d6) = ($d4, $d5, $d6)")
        # println("(e4, e5, e6) = ($e4, $e5, $e6)")
        # println("(p4, p5, p6) = ($p4, $p5, $p6)")
        # println("(t4, t5, t6) = ($t4, $t5, $t6)")
        # println("(q4, q5, q6) = ($q4, $q5, $q6)")
        # println("(r4, r5, r6) = ($r4, $r5, $r6)")

        for ip in get_integration_points(element)

            # shape functions

            xi, eta = ip.coords
            N1 = 2.0*(1.0-xi-eta)*(0.5-xi-eta)
            N2 = xi*(2.0*xi-1.0)
            N3 = eta*(2.0*eta-1.0)
            N4 = 4.0*xi*eta
            N5 = 4.0*eta*(1.0-xi-eta)
            N6 = 4.0*xi*(1.0-xi-eta)

            Hx = [
                1.5*(a6*N6 - a5*N5)
                b5*N5 + b6*N6
                N1 - c5*N5 - c6*N6
                1.5*(a4*N4 - a6*N6)
                b4*N4 + b6*N6
                N2 - c4*N4 - c6*N6
                1.5*(a5*N5 - a4*N4)
                b4*N4 + b5*N5
                N3 - c4*N4 - c5*N5]

            Hy = [
                1.5*(d6*N6 - d5*N5)
                -N1 + e5*N5 + e6*N6
                -b5*N5 - b6*N6
                1.5*(d4*N4 - d6*N6)
                -N2 + e4*N4 + e6*N6
                -b4*N4 - b6*N6
                1.5*(d5*N5 - d4*N4)
                -N3 + e4*N4 + e5*N5
                -b4*N4 - b5*N5]

            dHxdxi = [ p6*(1.0-2.0*xi) + (p5-p6)*eta
                       q6*(1.0-2.0*xi) - (q5+q6)*eta
                      -4.0 + 6.0*(xi+eta) + r6*(1.0-2.0*xi) - eta*(r5+r6)
                      -p6*(1.0-2.0*xi) + eta*(p4+p6)
                       q6*(1.0-2.0*xi) - eta*(q6-q4)
                      -2.0 + 6.0*xi + r6*(1.0-2.0*xi) + eta*(r4-r6)
                      -eta*(p5+p4)
                       eta*(q4-q5)
                      -eta*(r5-r4)]
            dHydxi = [ t6*(1.0-2.0*xi) + eta*(t5-t6)
                       1.0 + r6*(1.0-2.0*xi) - eta*(r5+r6)
                      -q6*(1.0-2.0*xi) + eta*(q5+q6)
                      -t6*(1.0-2.0*xi) + eta*(t4+t6)
                      -1.0 + r6*(1.0-2.0*xi) + eta*(r4-r6)
                      -q6*(1.0-2*xi) - eta*(q4-q6)
                      -eta*(t4+t5)
                       eta*(r4-r5)
                      -eta*(q4-q5)]
            dHxdeta = [-p5*(1.0-2.0*eta) - xi*(p6-p5)
                        q5*(1.0-2.0*eta) - xi*(q5+q6)
                        -4.0 + 6.0*(xi+eta) + r5*(1.0-2.0*eta) - xi*(r5+r6)
                        xi*(p4+p6)
                        xi*(q4-q6)
                        -xi*(r6-r4)
                        p5*(1.0-2.0*eta) - xi*(p4+p5)
                        q5*(1.0-2.0*eta) + xi*(q4-q5)
                        -2.0 + 6.0*eta + r5*(1.0-2.0*eta+xi*(r4-r5))]
            dHydeta = [-t5*(1.0-2.0*eta) - xi*(t6-t5)
                        1.0 + r5*(1.0-2.0*eta) - xi*(r5+r6)
                       -q5*(1.0-2.0*eta) + xi*(q5+q6)
                        xi*(t4+t6)
                        xi*(r4-r6)
                       -xi*(q4-q6)
                        t5*(1.0-2.0*eta) - xi*(t4+t5)
                       -1.0 + r5*(1.0-2.0*eta) + xi*(r4-r5)
                       -q5*(1.0-2.0*eta) - xi*(q4-q5)]

            # D11 = D22 = E*t^3/(12*(1-nu^2)) * 1.0
            # D21 = D22 = E*t^3/(12*(1-nu^2)) * nu
            # D66 = E*t^3/(12*(1-nu^2)) * 0.5*(1.0-nu)
            # Db = [D11 D12 0; D21 D22 0; 0 0 D66]

            # Material properties

            detJ = element(ip, time, Val{:detJ})
            E = element("youngs modulus", ip, time)
            t = element("thickness", ip, time)
            nu = element("poissons ratio", ip, time)

            Db = E*t^3/(12*(1-nu^2)) * [
                1.0  nu  0.0
                nu  1.0  0.0
                0.0 0.0  0.5*(1-nu)]

            # a11 = [ y3*p6            0                  -4*y3
            #        -y3*p6            0                   2*y3
            #         y3*p5           -y3*q5               y3*(2-r5)]
            # a12 = [-y3*p6            0                  -2*y3
            #         y3*p6            0                   4*y3
            #         y3*p4            y3*q4               y3*(r4-2)]
            # a13 = [ 0                0                   0
            #         0                0                   0
            #        -y3*(p4+p5)       y3*(q4-q5)          y3*(r4-r5)]
            # a21 = [-x2*t5            x23+x2*r5          -x2*q5
            #         0                x23                 0
            #         x23*t5           x23*(1-r5)          x23*q5
            # a22 = [  0               x3                  0
            #        x2*t4             x3+x2*r4           -x2*q4
            #       -x3*t4             x3*(1-r4)           x3*q4]
            # a23 = [x2*t5             x2*(r5-1)          -x2*q5
            #       -x2*t4             x2*(r4-1)          -x2*q4
            #       -x23*t5+x3*t4     -x23*r5-x3*r4-x2     x3*q4+x23*q5]
            # a31 = [-x3*p6-x2*p5      x2*q5+y3           -4*x23+x2*r5
            #        -x23*p6           y3                  2*x23
            #         x23*p5+y3*t5    -x23*q5+(1-r5)*y3   (2-r5)*x23+y3*q5]
            # a32 = [x3*p6            -y3                  2*x3
            #        x23*p6+x2*p4     -y3+x2*q4           -4*x3+x2*r4
            #       -x3*p4+y3*t4       (r4-1)*y3-x3*q4    (2-r4)*x3-y3*q4]
            # a33 = [x2*p5                       x2*q5                        (r5-2)*x2
            #       -x2*p4                       x2*q4                        (r4-2)*x2
            #       -x23*p5+x3*p4-(t4+t5)*y3    -x23*q5-x3*q4+(r4-r5)*y3     -x23*r5-x3*r4+4*x2+(q5-q4)*y3]

            # a = [a11 a12 a13
            #      a21 a22 a23
            #      a31 a32 a33]

            B = 1.0/(x31*y12 - x12*y31) * [
                 y31*dHxdxi' + y12*dHxdeta'
                -x31*dHydeta' - x12*dHydeta'
                -x31*dHxdxi' - x12*dHxdeta' + y31*dHydxi' + y12*dHydeta']

            # println("B = $B")
            # R = 1/24*[2 1 1; 1 2 1; 1 1 2]
            # z = UniformScaling(0.0)
            # F = 1/24*[
            #     D11*R D12*R z
            #     D21*R D22*R z
            #     z     z     D66*R]

            Ke += ip.weight * B'*Db*B * detJ

            N = element(ip, time)

            # if problem.properties.geometric_stiffness
            #     # Assemble Kg
            # end

            if haskey(element, "distributed load")
                p = element("distributed load", ip, time)
                fe[1:3:end] += ip.weight * p * N' * detJ
            end

        end # integration points


        gdofs = get_gdofs(problem, element)
        add!(assembly.K, gdofs, gdofs, Ke)
        add!(assembly.f, gdofs, fe)

    end # elements

    return nothing

end
