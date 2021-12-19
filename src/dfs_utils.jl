#=
This file contains structs and functions used the branch and bound
implementation.
=#

###
### Structs to hold upper and lower bounds on the treewidth of a graph.
###

"""A struct to hold the current upper bound of the branch and bound search."""
mutable struct UpperBound
    tw::UInt16            # The best treewidth found so far.
    order::Vector{UInt16} # The best elimination order found so far.
    lock::SpinLock        # A lock to prevent multiple threads updating the tw at once.
end

"""Use a heuristic to get an initial upper bound on the treewidth."""
function initialise_upper_bound(G::lg.AbstractGraph)
    # initial_ub_order, initial_ub_tw = min_fill(G; seed=seed)
    initial_ub_order, initial_ub_tw = collect(1:lg.nv(G)), 500
    initial_ub = UpperBound(initial_ub_tw, initial_ub_order, SpinLock())
    println("The initial treewidth is ", initial_ub_tw)
    initial_ub
end

"""A struct to store best intermediate treewidths. (see FPBB Yaun 2011)"""
struct LowerBounds
    tabel::Dict{BitVector, UInt16}
    lock::SpinLock

    LowerBounds() = new(Dict{BitVector, UInt16}(), SpinLock())
end


###
### A struct to maintain the search state of a depth first search through elimination orders.
###

"""Represents the search state of a single instance of the branch and bound search."""
mutable struct DFSState
    graph::Graph                       # The graph corresponding to the current state.
    N::UInt16                          # The size of the original graph.

    intermediate_graph_key::BitVector  # A bit vector indicating which vertices of the original graph
                                       # have been eliminated.

    curr_order::Vector{UInt16}         # The current elimination order being constructed.
    ub::UpperBound                     # A struct holding the best elimination order found so far.
    lbs::LowerBounds                   # A struct holding the lower bounds found so far (See Yaun 2011)

    branches::Vector{Vector{UInt16}}   # A vector of vetrices to search for each node depth

    pq::MMDQueue                       #

    # Some varibales to record where pruning happends (for introspection)
    lbs_prune_at_depth::Vector{UInt64}
    mmd_prune_at_depth::Vector{UInt64}
    ub_prune_at_depth::Vector{UInt64}

    # Variables to store various quantities related to the search.
    finish_time::Float64
    timed_out::Bool
    space_covered::Float64 # TODO: Try this as a BigFloat
    nodes_visited::UInt64
    depth::UInt16
    rng::MersenneTwister
end

"""DFSState constructor."""
function DFSState(g::Graph, ub::UpperBound, lbs::LowerBounds, finish_time::Float64, seed::Int=42)
    N = g.num_vertices
    rng = MersenneTwister(seed)
    initial_curr_order = collect(0x0001:N)

    branches = [Vector{UInt16}(undef, N-i) for i = 0:N-2]

    pq = MMDQueue(g.vertices, g.degree)

    DFSState(g, N, falses(N),
             initial_curr_order, ub, lbs, branches, pq,
             zeros(UInt64, N), zeros(UInt64, N), zeros(UInt64, N), 
             finish_time, false, 0.0, 0, 1, rng)
end

"""Return a copy of the given DFSState."""
function Base.copy(s::DFSState, seed::Int=42)
    DFSState(deepcopy(s.graph), s.N, copy(s.intermediate_graph_key),
             copy(s.curr_order), s.ub, s.lbs, deepcopy(s.branches), deepcopy(s.pq),
             zeros(UInt64, s.N), zeros(UInt64, s.N), zeros(UInt64, s.N),
             s.finish_time, false, 0.0, 0, s.depth, MersenneTwister(seed))
end

function initialise_dfsstates(G::lg.AbstractGraph, 
                                ub::UpperBound, 
                                finish_time::Float64, 
                                seed::Integer, 
                                n::Integer)
    # Create a search state for each thread.
    lbs = LowerBounds()
    state = DFSState(Graph(G), ub, lbs, finish_time, seed)
    # _, curr_tw = remove_simplicial_vertices!(state, 0)
    [copy(state, seed + t) for t = 1:n], 0x0000
end

"""Prints some info about the given DFSState."""
function Base.show(io::IO, s::DFSState)
    compact = get(io, :compact, false)
    if !compact
        println(io, "The graph size is {nodes, edges} = ", s.graph)
        println(io, "The best treewidth found was ", s.ub.tw)
        println(io, "The number of nodes visted was ", s.nodes_visited)
    else
        show(io, (s.ub.tw, s.ub.order))
    end
end


###
### A struct to hold the results of a depth first search.
###

"""Holds the results of a depth first search"""
struct DFSReport
    states::Vector{DFSState}    # A vector containing the states returned from each thread used.
    best_tw::Int                # The best treewidth found
    best_order::Vector{Int}     # The best elimination order found

    # Variables to indicate how much of the 
    # search space was covered.
    nodes_visited::UInt64
    space_covered::Float64

    # Time measurements.
    duration::Float64
    actual_time::Float64
    heuristic_time::Float64
    dfs_time::Float64

    # Vectors indicating where pruning happened.
    ub_pruned::Vector{UInt64}
    lbs_pruned::Vector{UInt64}
    mmd_pruned::Vector{UInt64}
end

"""Constructor for a DFSReport"""
function DFSReport(states::Vector{DFSState}, allocated_time, total_time, ub_time, dfs_time)
    best_tw = states[1].ub.tw
    best_order = copy(states[1].ub.order)

    total_nodes_visited = sum([state.nodes_visited for state in states])
    total_space_covered = sum([s.timed_out ? s.space_covered : 1/length(states) for s in states])

    ub_pruned = sum([state.ub_prune_at_depth for state in states])
    lbs_pruned = sum([state.lbs_prune_at_depth for state in states])
    mmd_pruned = sum([state.mmd_prune_at_depth for state in states])

    DFSReport(states, 
              best_tw, 
              best_order, 
              total_nodes_visited, 
              total_space_covered,
              allocated_time,
              total_time,
              ub_time,
              dfs_time,
              ub_pruned,
              lbs_pruned,
              mmd_pruned)
end

"""Prints the results of a branch-and-bound depth first search."""
function Base.show(io::IO, r::DFSReport)
    compact = get(io, :compact, false)
    if !compact
        println(io, "")
        println(io, "The graph size was {nodes, edges} = ", r.states[1].graph)
        println(io, "The best treewidth found was ", r.best_tw)
        println(io, "The number of nodes visted was ", r.nodes_visited)
        println(io, "The fraction of the search space covered was (to machine precision) ", r.space_covered)
        println(io, "Time:")
        println(io, " allocated = ", r.duration, " actual = ", r.actual_time)
        println(io, " heuristic = ", r.heuristic_time, " dfs = ", r.dfs_time)
    else
        show(io, (r.best_tw, r.best_order))
    end
end


###
### functions to alter the graph while searching through elimination orders
###

"""
eliminate!(state::DFSState, v::Int)

Eliminates vertex 'v' from the given state and update all its relevant variables.

Returns a vector of neighbours of v and the edges added when v was eliminated. These can
be used to restore the eliminated vertex.
"""
function Graphs.eliminate!(state::DFSState, v::UInt16)
    state.intermediate_graph_key[v] = true
    eliminate!(state.graph, v)
end

"""
    restore_last_eliminated!(state::DFSState,
                            edges_to_remove::Vector{Tuple{Int, Int}})

Restores a vertex 'v' which was eliminated to form the given 'state'.
"""
function Graphs.restore_last_eliminated!(state::DFSState,
                                   edges_to_remove::Vector{Tuple{UInt16, UInt16}})
    restore_last_eliminated!(state.graph, edges_to_remove)
    v = state.graph.vertices[state.graph.num_vertices]
    state.intermediate_graph_key[v] = false
    nothing
end

###
### 
###

"""Use the mmd heuristic on the given state graph"""
function mmd(state::DFSState)
    g = state.graph
    verts = vertices(g)

    initialise_mmdqueue!(state.pq, verts, g.degree)
    mmd(state.graph, state.pq)
end