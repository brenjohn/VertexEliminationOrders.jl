@testset "Heuristic tests.jl" begin
    # Load a graph to test on.
    g = graph_from_gr("../examples/test.gr")

    # Test the min fill heuristic.
    order, tw = min_fill(g)
    @test length(order) == nv(g)

    # Test the min width heuristic.
    order, tw = min_width(g)
    @test length(order) == nv(g)

    # Test sampling version of min width.
    order, tw = min_width_sampling(g, 10)
    @test length(order) == nv(g)

    # Test multi-threaded sampling version of min width.
    order, tw = min_width_mt_sampling(g, 10)
    @test length(order) == nv(g)
end