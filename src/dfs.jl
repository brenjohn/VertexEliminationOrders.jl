export myBB

###
### A struct to aid searching for elimination orders.
###

mutable struct UpperBound
    tw::Int
    order::Vector{Int}
    lock::SpinLock
end

struct LowerBounds
    tabel::Dict{Set{Int64}, Int64}
    lock::SpinLock

    LowerBounds() = new(Dict{Set{Int64}, Int64}(), SpinLock())
end

mutable struct DFSState
    graph::SimpleGraph{Int64}
    N::Int64
    labels::Vector{Int64}

    curr_order::Vector{Int64}
    ub::UpperBound
    lbs::LowerBounds

    lbs_prune_at_depth::Vector{UInt64}
    mmd_prune_at_depth::Vector{UInt64}
    ub_prune_at_depth::Vector{UInt64}

    finish_time::Float64
    space_covered::Float64 # TODO: Try this as a BigFloat
    nodes_visited::UInt64
    depth::Int
    rng::MersenneTwister
end


function DFSState(g::AbstractGraph, ub::UpperBound, lbs::LowerBounds, finish_time::Float64, seed::Int=42)
    N = nv(g)
    rng = MersenneTwister(seed)
    initial_curr_order = collect(1:N)

    DFSState(g, N, collect(1:N), 
             initial_curr_order, ub, lbs, 
             zeros(Int, N), zeros(Int, N), zeros(Int, N), 
             finish_time, 0.0 , 0, 1, rng)
end


function Base.show(io::IO, s::DFSState)
    compact = get(io, :compact, false)
    if !compact
        println(io, "The graph size was {nodes, edges} = ", s.graph)
        println(io, "The best treewidth found was ", s.ub.tw)
        println(io, "The number of nodes visted was ", s.nodes_visited)
    else
        show(io, (s.ub.tw, s.ub.order))
    end
end


struct DFSReport
    states::Vector{DFSState}
    best_tw::Int
    best_order::Vector{Int}

    nodes_visited::UInt64
    space_covered::Float64

    duration::Float64
    actual_time::Float64
    heuristic_time::Float64
    dfs_time::Float64

    ub_pruned::Vector{UInt64}
    lbs_pruned::Vector{UInt64}
    mmd_pruned::Vector{UInt64}
end


function DFSReport(states::Vector{DFSState}, allocated_time, total_time, ub_time, dfs_time)
    best_tw = states[1].ub.tw
    best_order = copy(states[1].ub.order)

    total_nodes_visited = sum([state.nodes_visited for state in states])
    total_space_covered = sum([state.space_covered for state in states])

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
### A branch and bound implementation
###

function myBB(G::AbstractGraph, duration::Float64, seed::Int=42)
    start_time = time()
    finish_time = start_time + duration
    # Use a heuristic to get an initial upper bound on the treewidth.
    initial_ub_order, initial_ub_tw = min_fill(G, seed)
    # initial_ub_order, initial_ub_tw = collect(1:nv(G)), nv(G)
    initial_ub = UpperBound(initial_ub_tw, initial_ub_order, SpinLock())
    heurisitc_finish_time = time()
    println("The initial treewidth is ", initial_ub_tw)

    # Create a search state for each thread.
    lbs = LowerBounds()
    states = [DFSState(deepcopy(G), initial_ub, lbs, finish_time, seed + t) for t = 1:nthreads()]

    # Randomly distribute the vertices of G amongst the threads.
    verts = shuffle(states[1].rng, collect(1:nv(G)))
    verts = [[verts[i] for i = j:nthreads():nv(G)] for j = 1:nthreads()]

    # Start each thread running a bb dfs.
    @threads for t = 1:nthreads()
        bb(states[t], verts[t], 0)
    end
    actual_finish_time = time()

    # Collect the results into a report.
    heuristic_time = heurisitc_finish_time - start_time
    actual_time = actual_finish_time - start_time
    dfs_time = actual_finish_time - heurisitc_finish_time
    DFSReport(states, duration, actual_time, heuristic_time, dfs_time)
end


function bb(state::DFSState, verts::Vector{Int}, curr_tw::Int)
    for (i, v) in enumerate(verts)
        if time() >= state.finish_time
            state.space_covered *= 1/length(verts)
            state.space_covered += (i-2)/length(verts) # -2 instead of -1 since i is incremented before check
            break
        end

        state.curr_order[state.depth] = state.labels[v]
        next_tw = max(curr_tw, degree(state.graph, v))

        v_neighbours, edges_added = modify_graph!(state, v)

        if curr_tw < state.ub.tw && check_lower_bounds!(state, curr_tw)
            dfs(state, next_tw)
        else
            if curr_tw >= state.ub.tw
                state.ub_prune_at_depth[state.depth] += 1
            end
        end

        unmodify_graph!(state, v, v_neighbours, edges_added)
    end
end


function dfs(state::DFSState, curr_tw::Int)
    state.depth += 1
    state.nodes_visited += 1

    if nv(state.graph) - 1 <= curr_tw
        if curr_tw < state.ub.tw
            lock(state.ub.lock)
            if curr_tw < state.ub.tw
                println("Thread $(threadid()): The treewidth is now ", curr_tw)
                state.curr_order[state.depth:end] = state.labels[1:nv(state.graph)]
                state.ub.order[:] = state.curr_order[:]
                state.ub.tw = curr_tw
            end
            unlock(state.ub.lock)
        end

    else
        verts = randperm(state.rng, nv(state.graph))
        I = sortperm([degree(state.graph, v) for v in verts])
        verts = verts[I]

        bb(state, verts, curr_tw)
    end

    state.depth -= 1
    state
end

function check_lower_bounds!(state::DFSState, curr_tw::Int)
    if mmd(state.graph) >= state.ub.tw
        state.mmd_prune_at_depth[state.depth] += 1
        return false
    end
    return true

    # state_key = Set(state.curr_order[1:state.depth])
    # state_lb = get(state.lbs.tabel, state_key, nothing)

    # if state_lb === nothing
    #     lock(state.lbs.lock)
    #     state_lb === nothing && (state.lbs.tabel[state_key] = curr_tw)
    #     unlock(state.lbs.lock)
    #     state_key = nothing
    #     return true

    # else
    #     if state_lb > curr_tw
    #         lock(state.lbs.lock)
    #         state_lb > curr_tw && (state.lbs.tabel[state_key] = curr_tw)
    #         unlock(state.lbs.lock)
    #         state_key = nothing
    #         return true
    #     else
    #         state.lbs_prune_at_depth[state.depth] += 1
    #         state_key = nothing
    #         return false
    #     end
    # end
end



###
### functions to alter the graph while searching through elimination orders
###

function modify_graph!(state::DFSState, v::Int64)
    G = state.graph
    n = nv(G)
    v_neighbours = copy(all_neighbors(G, v))

    # connect neighbours of v together
    edges_added = Tuple{Int, Int}[]
    for i = 1:length(v_neighbours)-1
        for j = i+1:length(v_neighbours)
            if add_edge!(G, v_neighbours[i], v_neighbours[j])
                append!(edges_added, [(v_neighbours[i], v_neighbours[j])])
            end
        end
    end

    # remove v from G
    rem_vertex!(G, v)

    # update the labels array to have the correct order
    tmp = state.labels[n]
    state.labels[n] = state.labels[v]
    state.labels[v] = tmp

    v_neighbours, edges_added
end


function unmodify_graph!(state::DFSState, v::Int64, v_neighbours::Vector{Int64}, edges_to_remove::Vector{Tuple{Int, Int}})
    G = state.graph
    add_vertex!(G)
    n = nv(G)

    n_neighbours = copy(all_neighbors(G, v))
    for u in n_neighbours
        rem_edge!(G, u, v)
        add_edge!(G, u, n)
    end

    for u in v_neighbours
        add_edge!(G, u, v)
    end

    for (ui, uj) in edges_to_remove
        rem_edge!(G, ui, uj)
    end

    # Swap the labels back
    tmp = state.labels[n]
    state.labels[n] = state.labels[v]
    state.labels[v] = tmp
end