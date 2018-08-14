# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMPlates.jl/blob/master/LICENSE

struct MindlinPlate <: FieldProblem end

function FEMBase.get_unknown_field_name(::Problem{MindlinPlate})
    return "displacement"
end

function FEMBase.assemble_elements!(problem::Problem{MindlinPlate},
                                    assembly::Assembly,
                                    elements::Vector{Element{Quad4}},
                                    time::Float64)

    for element in elements
        ndofs = 4*3
        Ke = zeros(ndofs, ndofs)
        fe = zeros(ndofs)
        for ip in get_integration_points(element)
            detJ = element(ip, time, Val{:detJ})
            E = element("youngs modulus", ip, time)
            t = element("thickness", ip, time)
            nu = element("poissons ratio", ip, time)
            gamma = 5/6
            if haskey(element, "shear correction factor")
                gamma = element("shear correction factor", ip, time)
            end
            # http://what-when-how.com/the-finite-element-method/fem-for-plates-and-shells-finite-element-method-part-1/
            G = E/(2.0 * (1.0 - nu))
            Db = E*t^3/(12*(1-nu^2)) * [
                1.0  nu  0.0
                nu  1.0  0.0
                0.0 0.0 1-nu]
            Ds = G*t*gamma*[1.0 0.0; 0.0 1.0]
            B_I = zeros(3, 12)
            N = element(ip, time)
            Grad = element(ip, time, Val{:Grad})
            for j=1:4
                B_I[1, (j-1)*3+3] = Grad[1,j]
                B_I[2, (j-1)*3+2] = Grad[2,j]
                B_I[3, (j-1)*3+2] = Grad[1,j]
                B_I[3, (j-1)*3+3] = Grad[2,j]
            end
            B_O = zeros(2, 12)
            for j=1:4
                B_O[1, (j-1)*3+1] = Grad[1,j]
                B_O[1, (j-1)*3+3] = N[j]
                B_O[2, (j-1)*3+1] = Grad[2,j]
                B_O[2, (j-1)*3+2] = N[j]
            end
            z = zeros(3, 2)
            D = [Db z; z' Ds]
            K1 = ip.weight * B_I'*Db*B_I * detJ
            K2 = ip.weight * B_O'*Ds*B_O * detJ
            Ke += (K1 + K2)
        end
        gdofs = get_gdofs(problem, element)
        add!(assembly.K, gdofs, gdofs, Ke)
        add!(assembly.f, gdofs, fe)
    end

end
