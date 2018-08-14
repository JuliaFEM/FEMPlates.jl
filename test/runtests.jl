# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMPlates.jl/blob/master/LICENSE

using Base.Test
using FEMPlates

@testset "FEMPlates.jl" begin
    @testset "test_mindlin_stiffness_matrix" begin
        include("test_mindlin_stiffness_matrix.jl")
    end
end
