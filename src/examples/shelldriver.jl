# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMPlates.jl/blob/master/LICENSE

using JuliaFEM
using FEMPlates

include("cylindrical_domain.jl")

function get_model()

    field_elements, bc_left, bc_right, bc_bottom, bc_top = get_rectangle()

    function load(element, ip, time)
        x,y = element("geometry", ip, time)
        return cos(2.0*y)
    end

    normals = Dict(1 => [0.0,0.0,1.0], 2 => [0.0,0.0,1.0], 3=> [0.0,0.0,1.0], 4 => [0.0,0.0,1.0])
    for element in field_elements
        for j in get_connectivity(element)
            normals[j] = [0.0,0.0,1.0]
        end
    end
    
    update!(field_elements, "youngs modulus", 1.0)
    update!(field_elements, "poissons ratio", 1/3)
    update!(field_elements, "thickness", 0.01)
    update!(field_elements, "radius of curvature", 1.0)
    update!(field_elements, "distributed load", load)
    update!(field_elements, "normal", normals)

    shell = Problem(Shell, "test problem", 5)
    shell_1elem = Problem(Shell, "test problem", 5)

    add_elements!(shell, field_elements)
    add_elements!(shell_1elem, [first(field_elements)])

    update!(bc_left, "displacement 1", 0.0)
    update!(bc_left, "displacement 4", 0.0)

    update!(bc_bottom, "displacement 2", 0.0)
    update!(bc_bottom, "displacement 5", 0.0)

    update!(bc_top, "displacement 1", 0.0)
    update!(bc_top, "displacement 3", 0.0)
    update!(bc_top, "displacement 4", 0.0)

    println("Assembly finished succesfully. Running static analysis, please wait.")

    shell.properties.debug = true

    return shell, bc_left, bc_right, bc_bottom, bc_top

end

function set_bc_1(bc_left, bc_right, bc_bottom, bc_top)

    update!(bc_right, "displacement 1", 0.0)
    update!(bc_right, "displacement 2", 0.0)
    update!(bc_right, "displacement 3", 0.0)
    update!(bc_right, "displacement 4", 0.0)
    update!(bc_right, "displacement 5", 0.0)

end

function set_bc_2(bc_left, bc_right, bc_bottom, bc_top)

end

function set_bc_3(bc_left, bc_right, bc_bottom, bc_top)

    update!(bc_right, "displacement 3", 0.0)

end

shell, bc_left, bc_right, bc_bottom, bc_top = get_model()
set_bc_1(bc_left, bc_right, bc_bottom, bc_top)
bc = Problem(Dirichlet, "fixed", 5, "displacement")
add_elements!(bc, bc_left)
add_elements!(bc, bc_right)
add_elements!(bc, bc_bottom)
add_elements!(bc, bc_top)

analysis = Analysis(Linear)
add_problems!(analysis, [shell, bc])
run!(analysis)

#for debugging:
#assemble!(shell_1elem, 0.0)
#K = full(shell_1elem.assembly.K)
#display(K)

X = shell("geometry", 0.0)
#u = shell.assembly.u
u = shell("displacement", 0.0)
N = length(X)
x = [X[i][1] for i=1:N]
y = [X[i][2] for i=1:N]
W = [u[i][3] for i=1:N]
d = []

Wref = [0,0]

Wref = [652.71419,652.13786,650.40969,647.53205,643.50882,638.34534,
632.04822, 624.62517, 616.08476, 606.43623, 595.68950, 583.85554,
570.94748, 556.98276, 541.98663, 525.99702, 509.07008, 491.28471,
472.74250, 453.55751, 433.82807, 413.58172, 392.68633, 370.72753,
346.86945, 319.74441, 287.46052, 247.87009, 199.29528, 141.93287,
80.111326, 25.381504, 0.0]

using Plots
surface(x, y, W)
gui()
