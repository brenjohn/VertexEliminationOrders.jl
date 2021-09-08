using VertexEliminationOrders
using Test

@testset "VertexEliminationOrders.jl" begin
    # Test the graph io.
    g = graph_from_gr("../examples/test.gr")
    @test nv(g) == 27
    @test ne(g) == 36
end

include("heuristic_tests.jl")
include("branch_and_bound_tests.jl")