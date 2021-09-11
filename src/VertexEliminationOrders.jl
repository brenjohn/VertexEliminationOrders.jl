module VertexEliminationOrders

using Base.Threads
using Random
using Reexport
@reexport using LightGraphs

include("utils.jl")
include("treewidth_heuristics.jl")
include("dfs_utils.jl")
include("dfs.jl")

end