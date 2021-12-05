@testset "Branch and bound tests.jl" begin
    # Load a graph to test on.
    g = graph_from_gr("../examples/test.gr")

    # Test the main branch and bound algorithm.
    dfs_report = branch_and_bound_dfs(g, 2.0)
    order = dfs_report.best_order
    @test length(order) == VertexEliminationOrders.lg.nv(g)
    @test dfs_report.best_tw == find_treewidth_from_order(g, order)


    
    # Create a DFSState to test removing and restoring simplicial vertices.
    g = graph_from_gr("../examples/sycamore_53_20_0.gr")
    state_i = VEO.DFSState(Graph(g), 
                       VEO.UpperBound(0, Int[], Base.Threads.SpinLock()), 
                       VEO.LowerBounds(), 0.0, 42)
    state_f = deepcopy(state_i)

    # # Remove and restore simplicial vertices.
    # removed_vertices, curr_tw = VEO.remove_simplicial_vertices!(state_f, 0)
    # VEO.restore_simplicial_vertices!(state_f, removed_vertices)

    # # Check if the graph is returned to its original state.
    # @test state_f.graph == state_i.graph
    # @test state_f.labels == state_i.labels
    # @test state_f.c_map == state_i.c_map
end