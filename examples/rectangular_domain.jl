# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMPlates.jl/blob/master/LICENSE

"""
    get_unit_square(nel_x=10, nel_y=10; lx=1.0, ly=1.0)
Return elements and boundary elements for a standard test problem, unit square.
"""
function get_rectangle(nel_x=32, nel_y=32; lx=1.0, ly=pi/4)

    nnodes_x = nel_x+1
    nnodes_y = nel_y+1
    nnode = nnodes_x*nnodes_y
    nodemap = reshape(1:nnode, nnodes_x, nnodes_y)
    nodes_1 = vec(nodemap[1:nnodes_x-1, 1:nnodes_y-1])
    nodes_2 = vec(nodemap[2:nnodes_x, 1:nnodes_y-1])
    nodes_3 = vec(nodemap[2:nnodes_x, 2:nnodes_y])
    nodes_4 = vec(nodemap[1:nnodes_x-1, 2:nnodes_y])

    X = Dict{Int64, Vector{Float64}}()

    # Create nodes
    nid = 1
    for y in linspace(0, ly, nnodes_y)
        for x in linspace(0, lx, nnodes_x)
            X[nid] = [x, y]
            nid += 1
        end
    end

    # Create elements for volume
    field_elements = []
    for c in zip(nodes_1, nodes_2, nodes_3, nodes_4)
        push!(field_elements, Element(Quad4, collect(c)))
    end

    # add boundary elements to the left side of domain
    bc_left = []
    bc_right = []
    bc_top = []
    bc_bottom = []
    nids = 1:nel_x+1:(nel_x+1)*(nel_y+1) # LEFT and RIGHT
    for c in zip(nids[1:end-1], nids[2:end])
        connectivity = collect(c)
        push!(bc_left, Element(Seg2, connectivity))
        push!(bc_right, Element(Seg2, connectivity+nel_x))
    end
    nids = 1:nel_x+1 # BOTTOM and TOP
    for c in zip(nids[1:end-1], nids[2:end])
        connectivity = collect(c)
        push!(bc_bottom, Element(Seg2, connectivity))
        push!(bc_top, Element(Seg2, connectivity+(nel_x+1)*nel_y))
    end

    update!(field_elements, "geometry", X)
    update!(bc_left, "geometry", X)
    update!(bc_right, "geometry", X)
    update!(bc_bottom, "geometry", X)
    update!(bc_top, "geometry", X)

    return field_elements, bc_left, bc_right, bc_bottom, bc_top
end
