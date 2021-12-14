export branch_and_bound_dfs

#=
This file contains an implementation of an anytime, branch and bound, depth first
search for computing the treewidth of a graph. If the search doesn't complete within
the allocated time, an upper bound is resturned representing the best solution the 
algorithm was able to find.

The algorithm is outlined by Gogate in the 2004 paper "A complete anytime algorithm
for treewidth".

Some of the methods mentioned in Yaun's 2011 paper "A Fast Parallel Branch and Bound 
Algorithm for Treewidth" and Bodlaender's "Treewidth Computations 1 & 2" papers are 
also used.

Some methods were ommitted but will hopefully be added eventually.
=#

#=
TODO: There seems to be a memory leak involving multiple threads and
storing the tabel of lowerbounds (when one of them is removed the memory
leak goes away). The leak results in the program OOM-ing when ran long enough.
=#

###
### A branch and bound implementation
###

"""
branch_and_bound_dfs(G::AbstractGraph, duration::Float63, seed::Int=42)
"""
function branch_and_bound_dfs(G::lg.AbstractGraph, duration::Float64, seed::Int=42)
    start_time = time()
    finish_time = start_time + duration
    
    # Use a heuristic to get an initial upper bound on the treewidth.
    # initial_ub_order, initial_ub_tw = min_fill(G; seed=seed)
    initial_ub_order, initial_ub_tw = collect(1:lg.nv(G)), 500
    initial_ub = UpperBound(initial_ub_tw, initial_ub_order, SpinLock())
    heurisitc_finish_time = time()
    println("The initial treewidth is ", initial_ub_tw)

    # Create a search state for each thread.
    lbs = LowerBounds()
    state = DFSState(Graph(G), initial_ub, lbs, finish_time, seed)
    # _, curr_tw = remove_simplicial_vertices!(state, 0)
    curr_tw = 0x0000
    states = [copy(state, seed + t) for t = 1:nthreads()]

    # Randomly distribute the vertices of G amongst the threads.
    verts = shuffle(states[1].rng, copy(state.graph.vertices))
    verts = [UInt16[verts[i] for i = j:nthreads():state.graph.num_vertices] for j = 1:nthreads()]

    # Start each thread running a bb-dfs.
    @threads for t = 1:nthreads()
        _bb(states[t], curr_tw, verts[t])
    end
    actual_finish_time = time()

    # Collect the results into a report.
    heuristic_time = heurisitc_finish_time - start_time
    actual_time = actual_finish_time - start_time
    dfs_time = actual_finish_time - heurisitc_finish_time
    DFSReport(states, duration, actual_time, heuristic_time, dfs_time)
end


function _bb(state::DFSState, curr_tw::UInt16, verts::Vector{UInt16})
    for v in verts
        bb(state, curr_tw, v)
    end
    nothing
end


"""Performs the branch and bound step"""
function bb(state::DFSState, curr_tw::UInt16, v::UInt16)
    # Add the next vertex to the current order and update the treewidth.
    state.curr_order[state.depth] = v
    next_tw = max(curr_tw, degree(state.graph, v))

    # Prune unless the current treewidth is less than the current upper bound.
    if next_tw < state.ub.tw #&& check_lower_bounds!(state, curr_tw)

        # Eliminate the vertex from the current graph.
        edges_added = eliminate!(state, v)
        state.nodes_visited += 1
        
        if check_lower_bounds!(state, next_tw)
            dfs(state, next_tw)
        end
        
        # Restore the graph to its original state before returning.
        restore_last_eliminated!(state, edges_added)
    else
        if curr_tw >= state.ub.tw
            state.ub_prune_at_depth[state.depth] += 1
        end
    end
end


"""Performs the recursive depth first search step"""
function dfs(state::DFSState, curr_tw::T) where T <: UInt16
    state.depth += 1

    # Recursive base case: an n-vertex graph can't have a treewidth 
    # larger than n-1.
    if state.graph.num_vertices - 0x0001 <= curr_tw
        if curr_tw < state.ub.tw
            lock(state.ub.lock)
            if curr_tw < state.ub.tw
                println("Thread $(threadid()): The treewidth is now ", curr_tw)
                state.curr_order[state.depth:end] .= vertices(state.graph)[:]
                copy!(state.ub.order, state.curr_order)
                state.ub.tw = curr_tw
            end
            unlock(state.ub.lock)
        end

    else
        # Reduce the graph if possible before moving to the next branch and bound step.
        # removed_vertices, curr_tw = remove_simplicial_vertices!(state, curr_tw)

        # Compute the vertices that need to be iterated over
        # in the next branch and bound step using theorem 6.1
        # in Gogate 2004. TODO: theorem 6.1 not implemented
        N = state.N - state.depth
        verts = state.branches[state.depth]
        copy!(verts, vertices(state.graph))

        # Loop over all candidates for the next vertex in the current order.
        for i = 0x0001:N
            v = verts[i]
            timed_out(state, N, i) && break

            bb(state, curr_tw, v)
        end

        # Restore the graph to its original state before returning.
        # restore_simplicial_vertices!(state, removed_vertices)
    end

    state.depth -= 1
    nothing
end


###
### functions to aid the branch and bound depth first search
###

"""Checks if pruning can be performed due to lower bounds"""
function check_lower_bounds!(state::DFSState, curr_tw::UInt16)
    # Compute a lower bound on the remaining graph and prune if
    # it is equal or larger than the current upper bound. 
    if mmd(state.graph) >= state.ub.tw
        state.mmd_prune_at_depth[state.depth] += 1
        return false
    end
    return true

    # # Below the pruning rule used by Yaun's tabel of lower bounds
    # # is used.
    # state_key = state.intermediate_graph_key
    # state_lb = get(state.lbs.tabel, state_key, nothing)

    # if state_lb === nothing
    #     lock(state.lbs.lock)

    #     if !haskey(state.lbs.tabel, state_key)
    #         state.lbs.tabel[state_key] = curr_tw
    #     elseif state.lbs.tabel[state_key] > curr_tw
    #         state.lbs.tabel[state_key] = curr_tw
    #     end

    #     unlock(state.lbs.lock)
    #     return true

    # else
    #     if state_lb > curr_tw
    #         lock(state.lbs.lock)
    #         state_lb > curr_tw && (state.lbs.tabel[state_key] = curr_tw)
    #         unlock(state.lbs.lock)
    #         return true
    #     else
    #         state.lbs_prune_at_depth[state.depth] += 1
    #         return false
    #     end
    # end
end


"""Implements the simplicial vertex rule for graph reduction (see Gogate 2004)"""
function remove_simplicial_vertices!(state::DFSState, curr_tw::Int)
    removed_vertices = Tuple{Int, Vector{Int}, Vector{Tuple{Int, Int}}}[]

    v = next_simplicial_vertex(state)
    while v !== nothing && nv(state.graph) - 1 > curr_tw
        curr_tw = max(curr_tw, degree(state.graph, v))
        state.curr_order[state.depth] = state.labels[v]
        state.depth += 1

        Nv, edges_added = eliminate!(state, v)
        push!(removed_vertices, (v::Int, Nv, edges_added))

        v = next_simplicial_vertex(state)
    end

    removed_vertices, curr_tw
end

"""Restores the given simplicial vertices that were removed from the graph."""
function restore_simplicial_vertices!(state::DFSState, 
                                      removed_vertices::Vector{Tuple{Int, Vector{Int}, Vector{Tuple{Int, Int}}}})
    state.depth -= length(removed_vertices)
    for (v, Nv, edges_added) in Iterators.reverse(removed_vertices)
        un_eliminate!(state, v, Nv, edges_added)
    end
end

"""Returns the index of a simplicial vertex which can be removed."""
function next_simplicial_vertex(state::DFSState)
    findfirst(c_v -> c_v == 0, state.c_map[1:nv(state.graph)])
end

"""Check if the algorithm has ran out of time"""
function timed_out(state::DFSState, N::UInt16, i::Integer)
    if time() >= state.finish_time
        if state.timed_out
            state.space_covered *= 1/N
            state.space_covered += (i-2)/N # -2 instead of -1 since i is incremented before check
        else
            state.timed_out = true
        end
        return true
    end
    false
end