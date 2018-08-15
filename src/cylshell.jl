# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMPlates.jl/blob/master/LICENSE

mutable struct CylShell <: FieldProblem
    debug :: Bool
end

function CylShell()
    return CylShell(false)
end

function FEMBase.get_unknown_field_name(::Problem{CylShell})
    return "displacement"
end

function FEMBase.assemble_elements!(problem::Problem{CylShell},
                                    assembly::Assembly,
                                    elements::Vector{Element{Quad4}},
                                    time::Float64)

    for element in elements
        ndofs = 4*5
        Ke = zeros(ndofs, ndofs)
        fe = zeros(ndofs)
        for ip in get_integration_points(element)
            detJ = element(ip, time, Val{:detJ})
            E = element("youngs modulus", ip, time)
            t = element("thickness", ip, time)
            nu = element("poissons ratio", ip, time)
            R = element("radius of curvature", ip, time)
            gamma = 5/6
            if haskey(element, "shear correction factor")
                gamma = element("shear correction factor", ip, time)
            end
            G = E/(2.0 * (1.0 + nu))
            Dm = E*t/(1-nu^2) * [
            1.0 nu 0.0
            nu 1.0 0.0
            0.0 0.0 (1-nu)/2
            ]
            Db = 0.0*E*t^3/(12*(1-nu^2)) * [
                1.0  nu  0.0
                nu  1.0  0.0
                0.0 0.0 (1-nu)/2]
            Ds = G*t*gamma*[1.0 0.0; 0.0 1.0]
            B_m = zeros(3, ndofs)
            B_b = zeros(3, ndofs)
            B_s = zeros(2, ndofs)
            N = element(ip, time)
            xir,etar = ip.coords
            Nrxi = element((0.0, etar), time)
            Nreta = element((xir, 0.0), time)
            Grad = element(ip, time, Val{:Grad})
            for j=1:4
                B_m[1, (j-1)*5+1] = Grad[1,j]
                B_m[2, (j-1)*5+2] = Grad[2,j]
                B_m[2, (j-1)*5+3] = N[j]/R
                B_m[3, (j-1)*5+1] = Grad[2,j]
                B_m[3, (j-1)*5+2] = Grad[1,j]
                B_b[1, (j-1)*5+4] = Grad[1,j]
                B_b[2, (j-1)*5+5] = Grad[2,j]
                B_b[3, (j-1)*5+4] = Grad[2,j]
                B_b[3, (j-1)*5+5] = Grad[1,j]
                B_b[3, (j-1)*5+2] = Grad[1,j]/R
                B_s[1, (j-1)*5+4] = N[j]
                B_s[1, (j-1)*5+3] = Grad[1,j]
                B_s[2, (j-1)*5+5] = N[j]
                B_s[2, (j-1)*5+3] = Grad[2,j]
                B_s[2, (j-1)*5+2] = -N[j]/R
            end
            Km = ip.weight * B_m'*Dm*B_m * detJ
            Ks = ip.weight * B_s'*Ds*B_s * detJ
            Kb = ip.weight * B_b'*Db*B_b * detJ
            Ke += (Km + Ks + Kb)
            if problem.properties.debug
                update!(element, "membrane stiffness", time => Km)
                update!(element, "shear stiffness", time => Ks)
                update!(element, "bending stiffness", time => Kb)
            end
            if haskey(element, "distributed load")
                p = element("distributed load", ip, time)
                fe[3:5:end-2] += ip.weight * p * N' * detJ
            end
        end
        gdofs = get_gdofs(problem, element)
        add!(assembly.K, gdofs, gdofs, Ke)
        add!(assembly.f, gdofs, fe)
    end

end
