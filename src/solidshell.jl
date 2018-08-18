# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMPlates.jl/blob/master/LICENSE

"""
    SolidShell - Solid shell formulation

# References

"""
struct SolidShell <: FieldProblem
end

function FEMBase.get_unknown_field_name(::Problem{SolidShell})
    return "displacement"
end

function FEMBase.assemble_elements!(problem::Problem{SolidShell},
                                    assembly::Assembly,
                                    elements::Vector{Element{Hex8}},
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
