export estimate_treewidth

#=
This file contains an implementation of an anytime, branch-and-bound depth-first
search for computing the treewidth of a graph. If the search doesn't complete within
the allocated time, an upper bound is resturned representing the best solution the 
algorithm was able to find.

The basic algorithm is outlined by Gogate in the 2004 paper "A complete anytime 
algorithm for treewidth".

One of the heuristics mentioned in Yaun's 2011 paper "A Fast Parallel Branch and 
Bound Algorithm for Treewidth" is also used.
    
A discussion of relevant heuristics can be found in Bodlaender's "Treewidth 
Computations 1 & 2" papers but not yet used.
=#

###
### Interface for the branch-and-bound depth-first search.
###

"""
    estimate_treewidth(G::AbstractGraph, duration::Float63; seed::Int=42)

Uses a branch-and-bound depth-first search to compute a vertex elimination order
for the graph `g` with minimal treewidth. If the search does not complete within
the allocated time, the best treewidth and elimination order is returned.

# Positional arguments
- `G::SimpleGraph`: A LightGraphs SimpleGraph representing the graph to estimate
                    the treewidth for.
- `run_time::AbstractFloat`: The amount of time to run the algorithm for.

# Keyword arguments
- `seed::Integer=-1`: The seed used by flow cutter. Most be a non negative integer,
                      otherwise the seed is set by flow cutter.
"""
function estimate_treewidth(
                        G::lg.SimpleGraph, 
                        run_time::AbstractFloat;
                        seed::Integer=42
                        )
    start_time = time()
    finish_time = start_time + run_time
    
    # Compute an initial upper bound on the treewidth.
    initial_ub = initialise_upper_bound(G, seed)
    heurisitc_finish_time = time()

    # Create a depth-first search state for each thread and run the algorithm.
    states = initialise_dfsstates(G, initial_ub, finish_time, seed, nthreads())
    _estimate_treewidth!(states)
    actual_finish_time = time()

    # Return the best result and a report containing metadata.
    report = DFSReport(states, run_time, start_time, heurisitc_finish_time, actual_finish_time)
    report.tw, report.order, report
end

"""Run the algorithm with the given vector states."""
function _estimate_treewidth!(states::Vector{DFSState})::Nothing
    # Randomly distribute the vertices of G amongst the threads.
    state = states[1]
    verts = shuffle(state.rng, copy(state.graph.vertices))
    verts = [UInt16[verts[i] for i = j:nthreads():state.graph.num_vertices] for j = 1:nthreads()]

    # Start each thread running a bb-dfs.
    @threads for t = 1:nthreads()
        _bb!(states[t], 0x0000, verts[t])
    end
    nothing
end

"""
Call the branch and bound step on each vertex in the given array.
(Only used by the interface function)
"""
function _bb!(state::DFSState, curr_tw::UInt16, verts::Vector{UInt16})
    for v in verts
        bb!(state, curr_tw, v)
    end
    nothing
end

###
### Implementation of a branch-and-bound depth-first search of elimination orders.
###

"""Performs the branch and bound step."""
function bb!(state::DFSState, curr_tw::UInt16, v::UInt16)
    # Add the next vertex to the current order and update the treewidth.
    state.curr_order[state.depth] = v
    next_tw = max(curr_tw, degree(state.graph, v))

    # Prune unless the current treewidth is less than the current upper bound.
    if next_tw < state.ub.tw

        # Eliminate the vertex from the current graph.
        eliminate!(state, v)
        state.nodes_visited += 1
        
        if check_lower_bounds!(state, next_tw)
            dfs!(state, next_tw)
        end
        
        # Restore the graph to its original state before returning.
        restore_last_eliminated!(state)
    else
        if curr_tw >= state.ub.tw
            state.ub_prune_at_depth[state.depth] += 1
        end
    end
end


"""Performs the recursive depth first search step."""
function dfs!(state::DFSState, curr_tw::T) where T <: UInt16
    state.depth += 1

    # Recursive base case: an n-vertex graph can't have a treewidth 
    # larger than n-1.
    if state.graph.num_vertices - 0x0001 <= curr_tw
        if curr_tw < state.ub.tw
            lock(state.ub.lock)
            if curr_tw < state.ub.tw
                @info "Thread $(threadid()): The treewidth is now $curr_tw"
                state.curr_order[state.depth:end] .= vertices(state.graph)[:]
                copy!(state.ub.order, state.curr_order)
                state.ub.tw = curr_tw
            end
            unlock(state.ub.lock)
        end

    else
        # Copy the vertices which may be eliminated next.
        # (This is needed as vertices(graph) may mutate after
        # subsequent eliminations.)
        N = num_vertices(state.graph)
        verts = state.branches[state.depth]
        copy!(verts, vertices(state.graph))

        # Explore eliminating each of the next candidates.
        for i = 0x0001:N
            timed_out(state, N, i) && break
            v = verts[i]
            bb!(state, curr_tw, v)
        end
    end

    state.depth -= 1
    nothing
end


###
### functions to aid the branch-and-bound depth-first search
###

"""Checks if pruning can be performed due to lower bounds"""
function check_lower_bounds!(state::DFSState, curr_tw::UInt16)
    # Compute a lower bound on the remaining graph and prune if
    # it is equal or larger than the current upper bound.
    if state.lb_ub >= state.ub.tw
        state.lb_ub = mmd(state)
        if state.lb_ub >= state.ub.tw
            state.mmd_prune_at_depth[state.depth] += 1
            return false
        end
    end

    # Below the pruning rule used by Yaun's tabel of lower bounds
    # is used.
    state_key = state.intermediate_graph_key
    tabel = state.lbs.tabel
    slock = state.lbs.lock
    state_lb = get(tabel, state_key, nothing)

    if state_lb === nothing
        lock(slock)
        if !haskey(tabel, state_key) || tabel[state_key] > curr_tw
            tabel[state_key] = curr_tw
        end
        unlock(slock)
        return true

    else
        if curr_tw < state_lb
            lock(slock)
            tabel[state_key] = curr_tw
            unlock(slock)
            return true
        else
            state.lbs_prune_at_depth[state.depth] += 1
            return false
        end
    end
end

"""Check if the algorithm has ran out of time"""
function timed_out(state::DFSState, N::UInt16, i::Integer)
    if time() >= state.finish_time
        if state.timed_out
            state.space_covered *= 1/N
            state.space_covered += (i-1)/N
        else
            state.timed_out = true
        end
        return true
    end
    false
end