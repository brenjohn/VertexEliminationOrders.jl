module VertexEliminationOrders

using Base.Threads
using Random

# Sub-modules to support dfs and heuristics.
include("mmdqueues.jl"); include("graphs.jl")
using .MMDQueues, .Graphs
export Graph

# Utility functions for loading and analysing graphs.
include("utils.jl")

include("treewidth_heuristics.jl")
include("dfs_utils.jl")
include("dfs.jl")

end