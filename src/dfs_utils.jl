#=
This file contains structs and functions used the branch and bound
implementation.
=#

###
### Structs to aid searching for elimination orders.
###

"""A struct to hold the current upper bound of the branch and bound search."""
mutable struct UpperBound
    tw::Int            # The best treewidth found so far.
    order::Vector{Int} # The best elimination order found so far.
    lock::SpinLock     # A lock to prevent multiple threads updating the tw at once.
end

"""A struct to store best intermediate treewidths. (see FPBB Yaun 2011)"""
struct LowerBounds
    tabel::Dict{BitVector, Int64}
    lock::SpinLock

    LowerBounds() = new(Dict{BitVector, Int64}(), SpinLock())
end

"""Represents the search state of a single instance of the branch and bound search."""
mutable struct DFSState
    graph::SimpleGraph{Int64}          # The graph corresponding to the current state.
    N::Int64                           # The size of the original graph.
    labels::Vector{Int64}              # A vector mapping graph vertices to their labels.
    c_map::Vector{Int64}               # A vector mapping graph vertices to their 'cliqueness' 
                                       # (num edges that need to be added to make it simplicial.)
    intermediate_graph_key::BitVector  # A bit vector indicating which vertices of the original graph
                                       # have been eliminated.

    curr_order::Vector{Int64}          # The current elimination order being constructed.
    ub::UpperBound                     # A struct holding the best elimination order found so far.
    lbs::LowerBounds                   # A struct holding the lower bounds found so far (See Yaun 2011)

    # Some varibales to record where pruning happends (for introspection)
    lbs_prune_at_depth::Vector{UInt64}
    mmd_prune_at_depth::Vector{UInt64}
    ub_prune_at_depth::Vector{UInt64}

    # Variables to store various quantities related to the search.
    finish_time::Float64
    timed_out::Bool
    space_covered::Float64 # TODO: Try this as a BigFloat
    nodes_visited::UInt64
    depth::Int
    rng::MersenneTwister
end

"""DFSState constructor."""
function DFSState(g::AbstractGraph, ub::UpperBound, lbs::LowerBounds, finish_time::Float64, seed::Int=42)
    N = nv(g)
    rng = MersenneTwister(seed)
    initial_curr_order = collect(1:N)
    c_map = [cliqueness(g, v) for v in vertices(g)]

    DFSState(g, N, collect(1:N), c_map, falses(N),
             initial_curr_order, ub, lbs, 
             zeros(Int, N), zeros(Int, N), zeros(Int, N), 
             finish_time, false, 0.0, 0, 1, rng)
end

"""Return a copy of the given DFSState."""
function Base.copy(s::DFSState, seed::Int=42)
    DFSState(deepcopy(s.graph), s.N, copy(s.labels), copy(s.c_map), copy(s.intermediate_graph_key),
             copy(s.curr_order), s.ub, s.lbs,
             zeros(Int, s.N), zeros(Int, s.N), zeros(Int, s.N),
             s.finish_time, false, 0.0, 0, s.depth, MersenneTwister(seed))
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
function eliminate!(state::DFSState, v::Int64)
    G = state.graph
    n = nv(G)
    v_neighbours = copy(all_neighbors(G, v))

    # connect neighbours of v together
    edges_added = Tuple{Int, Int}[]
    for i = 1:length(v_neighbours)-1
        for j = i+1:length(v_neighbours)
            if add_edge!(G, v_neighbours[i], v_neighbours[j])
                @inbounds push!(edges_added, (v_neighbours[i], v_neighbours[j]))

                # Common neighbours of vi and vj have one less edge to add
                # when being eliminated after vi and ui are connected.
                @inbounds vi = v_neighbours[i]
                @inbounds vj = v_neighbours[j]
                Nvi = all_neighbors(G, vi)::Array{Int64, 1}
                Nvj = all_neighbors(G, vj)::Array{Int64, 1}
                for n in Nvi
                    if n in Nvj
                        @inbounds state.c_map[n] -= 1
                    end
                end

                # ui and vi are now neighbours so their cliqueness may increase.
                for n in Nvi
                    if !(n == vj) && !(has_edge(G, n, vj)::Bool)
                        @inbounds state.c_map[vi] += 1
                    end
                end
                for n in Nvj
                    if !(n == vi) && !(has_edge(G, n, vi)::Bool)
                        @inbounds state.c_map[vj] += 1
                    end
                end
            end
        end
    end

    # Removing v from G means it's also removed from its neighbour's neighbourhood, so their 
    # cliqueness may be reduced.
    for n in v_neighbours
        Nₙ = all_neighbors(G, n)::Array{Int64, 1}
        for u in Nₙ
            if !(u == v)
                if !has_edge(G, v, u)
                    @inbounds state.c_map[n] -= 1
                end
            end
        end
    end

    # remove v from G
    rem_vertex!(G, v)

    # update the labels array to have the correct order
    @inbounds begin
        tmp = state.labels[n]
        state.labels[n] = state.labels[v]
        state.labels[v] = tmp
        state.intermediate_graph_key[v] = true

        tmp = state.c_map[n]
        state.c_map[n] = state.c_map[v]
        state.c_map[v] = tmp
    end

    v_neighbours, edges_added
end

"""
    un_eliminate!(state::DFSState, 
                v::Int64, 
                v_neighbours::Vector{Int64}, 
                edges_to_remove::Vector{Tuple{Int, Int}})

Restores a vertex 'v' which was eliminated to form the given 'state'.
"""
function un_eliminate!(state::DFSState, 
                       v::Int64, 
                       v_neighbours::Vector{Int64}, 
                       edges_to_remove::Vector{Tuple{Int, Int}})
    G = state.graph
    add_vertex!(G)
    n = nv(G)

    @inbounds begin
        tmp = state.c_map[n]
        state.c_map[n] = state.c_map[v]
        state.c_map[v] = tmp
    end

    n_neighbours = copy(all_neighbors(G, v))
    for u in n_neighbours
        rem_edge!(G, u, v)
        add_edge!(G, u, n)
    end

    for u in v_neighbours
        add_edge!(G, u, v)
    end

    # Restoring v in G means it's neighbour's neighbourhoods increase, so their 
    # cliqueness may be increased.also removed from its 
    for n in v_neighbours
        Nₙ = all_neighbors(G, n)::Array{Int64, 1}
        for u in Nₙ
            if !(u == v)
                if !has_edge(G, v, u)
                    @inbounds state.c_map[n] += 1
                end
            end
        end
    end

    for (vi, vj) in edges_to_remove
        rem_edge!(G, vi, vj)

        # Common neighbours of vi and vj have one more edge to add
        # when being eliminated after vi and ui are connected.
        Nvi = all_neighbors(G, vi)::Array{Int64, 1}
        Nvj = all_neighbors(G, vj)::Array{Int64, 1}
        for n in Nvi
            if n in Nvj
                @inbounds state.c_map[n] += 1
            end
        end

        # vi and vi are no longer neighbours so their cliqueness may decrease.
        for n in Nvi
            if !(has_edge(G, n, vj)::Bool)
                @inbounds state.c_map[vi] -= 1
            end
        end
        for n in Nvj
            if !(has_edge(G, n, vi)::Bool)
                @inbounds state.c_map[vj] -= 1
            end
        end
    end

    # Swap the labels back
    @inbounds begin
        tmp = state.labels[n]
        state.labels[n] = state.labels[v]
        state.labels[v] = tmp
        state.intermediate_graph_key[v] = false
    end
    
    nothing
end