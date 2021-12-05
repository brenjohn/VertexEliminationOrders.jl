@testset "Heuristic tests" begin
    # Load a graph to test on.
    g = graph_from_gr("../examples/test.gr")

    # Test the min fill heuristic.
    order, tw = min_fill(g)
    @test length(order) == lg.nv(g)

    # Test min fill sampling
    order, tw = min_fill(g, 10)
    @test length(order) == lg.nv(g)

    # # Test the min width heuristic.
    # order, tw = min_width(g)
    # @test length(order) == lg.nv(g)

    # # Test sampling version of min width.
    # order, tw = min_width(g, 10)
    # @test length(order) == lg.nv(g)

    # Test maximum minimum degree heuristic for a lower bound.
    lb = mmd(g)
    @test lb <= tw
end