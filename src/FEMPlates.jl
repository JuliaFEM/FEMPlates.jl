# This file is a part of JuliaFEM.
# License is MIT: see https://github.com/JuliaFEM/FEMPlates.jl/blob/master/LICENSE

module FEMPlates

using Reexport
@reexport using FEMBase

include("mindlin.jl")
export MindlinPlate

include("dkt.jl")
export DKT

include("shell.jl")
export Shell

end
