export min_fill, min_width
export min_width_sampling
export min_width_mt_sampling
export mmd

#=
This file contains heuristic functions for computing upper and lower bounds
on the treewidth of a graph.
=#

###
### Min fill heuristic for upper bounds.
###

"""
    min_fill(g::AbstractGraph)

Use the min fill heuristic to find a vertex elimination order for g with
a minimal treewidth.
"""
function min_fill(g::lg.AbstractGraph; seed::Int=42)
    labels = collect(1:lg.nv(g))
    cliqueness_map = [cliqueness(g, v) for v in lg.vertices(g)]
    order = Array{Int, 1}(undef, lg.nv(g))
    tw = min_fill(MersenneTwister(seed), g, labels, cliqueness_map, order)
    order, tw
end

function min_fill(rng::AbstractRNG, 
                  g::lg.AbstractGraph, 
                  labels::Array{Int, 1}, 
                  cliqueness_map::Array{Int, 1}, 
                  order::Array{Int, 1})
    g = deepcopy(g)
    tw = 0

    i = 1
    n = lg.nv(g)
    while n > 0
        c_map = @view cliqueness_map[1:n]
        candidates = findall(isequal(minimum(c_map)), c_map)
        v = rand(rng, candidates)
        tw = max(tw, lg.degree(g, v))
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
function min_fill(g::lg.AbstractGraph, samples::Integer; seed::Integer=42)
    # Allocate memory for each thread.
    N = lg.nv(g)
    rngs = [MersenneTwister(seed + threadid()) for i = 1:nthreads()]
    orders = [Array{Int, 1}(undef, N) for _ = 1:nthreads()]
    graphs = [deepcopy(g) for i = 1:nthreads()]
    tws = Array{Int, 1}(undef, nthreads())

    cliqueness_map = [Array{Int, 1}(undef, N) for _ = 1:nthreads()]
    labels = [Array{Int, 1}(undef, N) for _ = 1:nthreads()]

    best_orders = [Array{Int, 1}(undef, N) for _ = 1:nthreads()]
    best_tws = fill(N, nthreads())


    @threads for i = 1:samples
        # Re-initialise cliqueness_map and label caches before next min fill call.
        @inbounds for v = 1:N
            cliqueness_map[threadid()][v] = cliqueness(g, v)
            labels[threadid()][v] = v
        end

        # call the min fill heuristic.
        tws[threadid()] = min_fill(rngs[threadid()], 
                                    graphs[threadid()],
                                    labels[threadid()],
                                    cliqueness_map[threadid()],
                                    orders[threadid()])

        # Update the best treewidth and order if a better one was found.
        if tws[threadid()] < best_tws[threadid()]
            best_orders[threadid()][:] = orders[threadid()][:]
            best_tws[threadid()] = tws[threadid()]
        end
    end

    # Return the best result.
    best = argmin(best_tws)
    best_orders[best], best_tws[best]
end


###
### Min Width Heuristic
###

# """
#     min_width(g::AbstractGraph, seed::Union{Nothing, Int})

# Use the min width heuristic to find a vertex elimination order for g with
# a minimal treewidth.
# """
# function min_width(g::AbstractGraph; seed::Union{Nothing, Int}=nothing)
#     labels = collect(1:nv(g))
#     order = Array{Int, 1}(undef, nv(g))
#     rng = seed isa Int ? MersenneTwister(seed) : MersenneTwister()
#     tw = min_width(rng, g, labels, order)
#     order, tw
# end

# function min_width(rng::AbstractRNG, 
#                    g::AbstractGraph, 
#                    labels::Array{Int, 1}, 
#                    order::Array{Int, 1})
#     g = deepcopy(g)
#     tw = 0

#     i = 1
#     n = nv(g)
#     while i <= n
#         d_map = degree(g)
#         candidates = findall(isequal(minimum(d_map)), d_map)
#         v = rand(rng, candidates)
#         tw = max(tw, degree(g, v))
#         order[i] = labels[v]
#         eliminate!(g, labels, v)
#         i += 1
#     end

#     tw
# end

# """
#     min_width(g::AbstractGraph, samples::Integer, seed::Integer=42)

# Samples a specified number of vertex elimination orders for g using the min fill 
# heuristic and returns the best one.

# Sampling is multi-threaded.
# """
# function min_width(g::AbstractGraph, samples::Int; seed::Int=42)
#     # Allocate memory for each thread.
#     orders = Array{Array{Int, 1}, 1}(undef, nthreads())
#     tws = Array{Int, 1}(undef, nthreads())

#     # Get each thread to sample elimination orders.
#     samples_per_thread = samples ÷ nthreads()
#     @threads for t = 1:nthreads()
#         orders[threadid()], tws[threadid()] = min_width_sampling(g, 
#                                                                 samples_per_thread + ((samples % nthreads()) > 0); 
#                                                                 seed = seed + t)
#     end

#     # return the best result.
#     best = argmin(tws)
#     orders[best], tws[best]
# end

# function min_width_sampling(g::AbstractGraph, N_samples::Int; seed::Int=42)
#     # Create an rng and allocate memory.
#     rng = MersenneTwister(seed)

#     best_order = Array{Int, 1}(undef, nv(g))
#     best_tw = nv(g)
    
#     curr_order = Array{Int, 1}(undef, nv(g))

#     # Sample the given number of elimination orders.
#     for i = 1:N_samples
#         curr_tw = min_width(rng, deepcopy(g), collect(1:nv(g)), curr_order)
#         if curr_tw < best_tw
#             best_tw = curr_tw
#             best_order[:] = curr_order[:]
#         end
#     end

#     # Return the best result.
#     best_order, best_tw
# end


###
### Lower bound heuristics (see Bodlaender 2011)
###

mmd(g::lg.SimpleGraph) = mmd(Graph(g))

"""Maximum Minimum Degree"""
function mmd(g::Graph)
    pq = MMDQueue(g.vertices, g.degree)
    mmd(g, pq)
end

function mmd(g::Graph, pq::MMDQueue)
    N = g.num_vertices - 0x0001 # max degree possible. 

    maxmin = 0x0000
    i = 1
    while maxmin < N
        @inbounds v_d = pull!(pq)
        maxmin = max(maxmin, v_d.second)

        for u in neighbours(g, v_d.first)
            decrement!(pq, u)
        end

        i += 1
        N -= 0x0001
    end

    maxmin
end

# """
# Maximum Minimum Degree plus

# This implementation is correct but slow
# """
# function mmd_plus(G::SimpleGraph{Int})
#     H = deepcopy(G)
#     maxmin = 0

#     while nv(H) > 2
#         v = argmin(degree(H))
#         maxmin = max(maxmin, degree(H, v))
        
#         ui = argmin([degree(H, n) for n in all_neighbors(H, v)])
#         u = all_neighbors(H, v)[ui]
#         merge_vertices!(H, [v, u])
#     end

#     maxmin
# end