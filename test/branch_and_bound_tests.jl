@testset "Branch and bound tests.jl" begin
    # Load a graph to test on.
    g = graph_from_gr("../examples/test.gr")

    # Test the main branch and bound algorithm.
    tw, order, report = estimate_treewidth(g, 2.0)
    @test length(order) == VertexEliminationOrders.lg.nv(g)
    @test tw == find_treewidth_from_order(g, order)
end