using VertexEliminationOrders
using Test

import LightGraphs as lg

import VertexEliminationOrders as VEO

@testset "VertexEliminationOrders.jl" begin
    # Test the graph io.
    g = graph_from_gr("../examples/test.gr")
    @test lg.nv(g) == 27
    @test lg.ne(g) == 36
end

include("graph_tests.jl")
include("heuristic_tests.jl")
include("branch_and_bound_tests.jl")