export min_fill, min_width
export min_width_sampling
export min_width_mt_sampling
export mmd, mmd_plus, mmd_plus_correct

###
### Min Fill Heuristic
###

"""
    min_fill(g::AbstractGraph)

Use the min fill heuristic to find a vertex elimination order for g.
"""
function min_fill(g::AbstractGraph, seed::Int=42)
    labels = collect(1:nv(g))
    cliqueness_map = [cliqueness(g, v) for v in vertices(g)]
    order = Array{Int, 1}(undef, nv(g))
    tw = min_fill(MersenneTwister(seed), g, labels, cliqueness_map, order)
    order, tw
end

function min_fill(rng::AbstractRNG, 
                  g::AbstractGraph, 
                  labels::Array{Int, 1}, 
                  cliqueness_map::Array{Int, 1}, 
                  order::Array{Int, 1})
    g = deepcopy(g)
    tw = 0

    i = 1
    n = nv(g)
    while n > 0
        c_map = @view cliqueness_map[1:n]
        candidates = findall(isequal(minimum(c_map)), c_map)
        v = rand(rng, candidates)
        tw = max(tw, degree(g, v))
        order[i] = labels[v]
        eliminate!(g, labels, cliqueness_map, v)
        n -= 1
        i += 1
    end

    tw
end


"""
    min_fill(g::AbstractGraph, samples::Integer, seed::Integer=42)

Samples a specified number of vertex elimination orders for g using the min fill 
heuristic and returns the best one.

Sampling is multi-threaded.
"""
function min_fill(g::AbstractGraph, samples::Integer, seed::Integer=42)
    N = nv(g)
    rngs = [MersenneTwister(seed + threadid()) for i = 1:nthreads()]
    orders = [Array{Int, 1}(undef, N) for _ = 1:nthreads()]
    graphs = [deepcopy(g) for i = 1:nthreads()]
    tws = Array{Int, 1}(undef, nthreads())

    cliqueness_map = [Array{Int, 1}(undef, N) for _ = 1:nthreads()]
    labels = [Array{Int, 1}(undef, N) for _ = 1:nthreads()]

    best_orders = Array{Array{Int, 1}}(undef, nthreads())
    best_tws = fill(nv(g), nthreads())


    @threads for i = 1:samples
        # TODO: This is slower for some reason, I think id and tw are being allocated on 
        # the heap and not the stack and so are shared between threads?:
        # id = Threads.threadid()::Int64
        # orders[id], tw = min_fill(rngs[id], graphs[id])
        # if tw < best_tws[id]
        #     best_orders[id] = orders[id]
        #     best_tws[id] = tw
        # end

        # Re-initialise cliqueness_map and label caches before next min fill call.
        @inbounds for v = 1:N
            cliqueness_map[threadid()][v] = cliqueness(g, v)
            labels[threadid()][v] = v
        end


        tws[threadid()] = min_fill(rngs[threadid()], 
                                    graphs[threadid()],
                                    labels[threadid()],
                                    cliqueness_map[threadid()],
                                    orders[threadid()])

        if tws[threadid()] < best_tws[threadid()]
            best_orders[threadid()] = orders[threadid()]
            best_tws[threadid()] = tws[threadid()]
        end
    end

    best = argmax(best_tws)
    best_orders[best], best_tws[best]
end


"""
    eliminate!(graph, labels, vertex)

Connects the neighbours of the given vertex into a clique before removing 
it from the graph.

The array of vertex labels are also updated to reflect the 
reording of vertex indices when the graph is updated.
"""
function eliminate!(g::AbstractGraph, labels, v)
    Nᵥ = all_neighbors(g, v)::Array{Int64, 1}
    for i = 1:length(Nᵥ)-1
        vi = Nᵥ[i]
        for j = i+1:length(Nᵥ)
            vj = Nᵥ[j]
            add_edge!(g, vi, vj)
        end
    end

    rem_vertex!(g, v)
    labels[v] = labels[end]
    pop!(labels)
    g
end

"""
    eliminate!(graph, labels, vertex)

Connects the neighbours of the given vertex into a clique before removing 
it from the graph.

The arrays of vertex labels and cliqueness are also updated to reflect the 
reording of vertex indices when the graph is updated.
"""
function eliminate!(g::AbstractGraph, labels, c_map, v)
    Nᵥ = all_neighbors(g, v)::Array{Int64, 1}
    for i = 1:length(Nᵥ)-1
        vi = Nᵥ[i]
        for j = i+1:length(Nᵥ)
            vj = Nᵥ[j]

            # Try add an edge connecting vi and ui. If successful, update `c_map`.
            edge_added = add_edge!(g, vi, vj)
            if edge_added
                # Common neighbours of vi and vj have one less edge to add
                # when being eliminated after vi and ui are connected.
                Nvi = all_neighbors(g, vi)::Array{Int64, 1}
                Nvj = all_neighbors(g, vj)::Array{Int64, 1}
                for n in Nvi
                    if n in Nvj
                        c_map[n] -= 1
                    end
                end

                # ui and vi are now neighbours so their cliqueness may increase.
                for n in Nvi
                    if !(n == vj) && !(has_edge(g, n, vj)::Bool)
                        c_map[vi] += 1
                    end
                end
                for n in Nvj
                    if !(n == vi) && !(has_edge(g, n, vi)::Bool)
                        c_map[vj] += 1
                    end
                end
            end
        end
    end

    # Removing v from G means it's also removed from its neighbour's neighbourhood, so their 
    # cliqueness may be reduced.
    for n in Nᵥ
        Nₙ = all_neighbors(g, n)::Array{Int64, 1}
        for u in Nₙ
            if !(u == v)
                if !has_edge(g, v, u)
                    c_map[n] -= 1
                end
            end
        end
    end

    rem_vertex!(g, v)
    c_map[v] = c_map[end]
    labels[v] = labels[end]
    # pop!(c_map)
    pop!(labels)
    g
end


"""
    cliqueness(G::AbstractGraph, v::Integer)

Return the number of edges that need to be added to `G` in order to make the neighborhood of 
vertex `v` a clique.
"""
function cliqueness(g::AbstractGraph, v::Integer)::Int
    neighborhood = all_neighbors(g, v)::Array{Int64, 1}
    count = 0
    for i in 1:length(neighborhood)-1
        for j in i+1:length(neighborhood)
            vi = neighborhood[i]
            ui = neighborhood[j]
            if !has_edge(g, ui, vi)::Bool
                count += 1
            end
        end
    end
    count
end

###
### Min Width Heuristic
###

function min_width(g::AbstractGraph, seed=nothing)
    labels = collect(1:nv(g))
    order = Array{Int, 1}(undef, nv(g))
    rng = seed isa Int ? MersenneTwister(seed) : MersenneTwister()
    tw = min_width(rng, g, labels, order)
    order, tw
end

function min_width(rng::AbstractRNG, 
                   g::AbstractGraph, 
                   labels::Array{Int, 1}, 
                   order::Array{Int, 1})
    g = deepcopy(g)
    tw = 0

    i = 1
    n = nv(g)
    while i <= n
        d_map = degree(g)
        candidates = findall(isequal(minimum(d_map)), d_map)
        v = rand(rng, candidates)
        tw = max(tw, degree(g, v))
        order[i] = labels[v]
        eliminate!(g, labels, v)
        i += 1
    end

    tw
end

function min_width_sampling(g::AbstractGraph, N_samples::Int)
    rng = MersenneTwister()

    best_order = Array{Int, 1}(undef, nv(g))
    best_tw = nv(g)
    
    curr_order = Array{Int, 1}(undef, nv(g))

    for i = 1:N_samples
        curr_tw = min_width(rng, deepcopy(g), collect(1:nv(g)), curr_order)
        if curr_tw < best_tw
            best_tw = curr_tw
            best_order[:] = curr_order[:]
        end
    end

    best_order, best_tw
end

function min_width_mt_sampling(g::AbstractGraph, samples_per_thread::Int)
    orders = Array{Array{Int, 1}, 1}(undef, nthreads())
    tws = Array{Int, 1}(undef, nthreads())

    @threads for _ = 1:nthreads()
        orders[threadid()], tws[threadid()] = min_width_sampling(g, samples_per_thread)
    end

    best = argmin(tws)

    orders[best], tws[best]
end


###
### Lower Bound heuristics
###


"""Maximum Minimum Degree"""
function mmd(G::SimpleGraph{Int})
    N = nv(G)
    degree_map = degree(G)
    maxmin = 0
    while maxmin < N # N > 2
        dv, v = findmin(degree_map)
        maxmin = max(maxmin, dv)

        @inbounds degree_map[v] = typemax(Int)
        for u in all_neighbors(G, v)
            @inbounds degree_map[u] -= 1
        end

        N -= 1
    end

    maxmin
end

"""
Maximum Minimum Degree +

TODO: This implementation is incorrect.
problem with neighbours arrays after vertices merged.
"""
function mmd_plus(G::SimpleGraph{Int})
    N = nv(G)
    degree_map = degree(G)
    maxmin = 0
    while maxmin < N # N > 2
        dv, v = findmin(degree_map)
        maxmin = max(maxmin, dv)

        w = 0
        dw = N
        for wi in all_neighbors(G, v) # not correct after vertices 'merged'
            @inbounds if degree_map[wi] < dw
                @inbounds dw = degree_map[wi]
                w = wi
            end
        end

        if w != 0 # This check can be removed if we assume the grraph is connected.
            for u in all_neighbors(G, v)
                if u in all_neighbors(G, w) 
                    @inbounds degree_map[u] -= 1
                else
                    @inbounds degree_map[w] += 1
                end
            end
        end
        degree_map[v] = typemax(Int)

        N -= 1
    end

    maxmin
end

"""
Maximum Minimum Degree plus

This implementation is correct but slow
"""
function mmd_plus_correct(G::SimpleGraph{Int})
    H = deepcopy(G)
    maxmin = 0

    while nv(H) > 2
        v = argmin(degree(H))
        maxmin = max(maxmin, degree(H, v))
        
        ui = argmin([degree(H, n) for n in all_neighbors(H, v)])
        u = all_neighbors(H, v)[ui]
        merge_vertices!(H, [v, u])
    end

    maxmin
end