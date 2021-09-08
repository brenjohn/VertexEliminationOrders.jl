@testset "Branch and bound tests.jl" begin
    # Load a graph to test on.
    g = graph_from_gr("../examples/test.gr")

    # Test the min fill heuristic.
    dfs_report = myBB(g, 2.0)

    order = dfs_report.best_order
    @test length(order) == nv(g)
end