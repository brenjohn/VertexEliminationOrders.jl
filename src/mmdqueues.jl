module MMDQueues

# TODO: I don't think the extra slot at the end of the queue is needed anymore.

#=
This module defines a custom priority queue optimised for implementing 
the mmd heuristic.

Retrieving the highest priority element and decrementing the priority
of elements by 1 are both order 1 operations in complextity. mmd does
not need an insertion operation so it can be omitted.

The queue is represented by a sorted array of (vertex, degree) pairs.
An array of keys or indices for the vertices is maintained for order
1 look up of a vertex in the queue.

A bag structure is also maintained to help keep the queue sorted
when the degree of a vertex is decremented. ie. the queue is 
partitioned into bags, with all vertices in bag d having degree d.
This is maintained by an array containing the positions of the left
most boundary of a bag. The vertices of degree d begin at bags[d+1]
inclusively and end at bags[d+2] exclusively.


            Bag 0 <----- empty Bags ------>  Bag 3
            |___|                            |___|
              |--Bag 1-|-------Bag 2-----------|--Bag 4-|--Bag 5-|
              |        |                       |        |        |       
              |        |                       |        |        |
queue: [ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ][ ]
        |  |  |
        |  |   \__
removed----       \--min
=#

export MMDQueue, pop!, decrement!

"""
A custom priority queue optimised for implementing the mmd heuristic.

Retrieving the highest priority element and decrementing the priority
of elements by 1 are both order 1 operations in complextity. mmd does
not need an insertion operation so it can be omitted.

See comments above docstring for implementation details.
"""
mutable struct MMDQueue
    queue::Vector{Pair{UInt16, UInt16}}
    keys::Vector{UInt16}
    bags::Vector{UInt16}
    min::UInt16
end

function MMDQueue(vertices::Vector{UInt16}, degrees::Vector{UInt16})
    I = sortperm(vertices; by=v -> degrees[v])
    N = length(vertices)

    queue = Array{Pair{UInt16, UInt16}, 1}(undef, N+1)
    keys = similar(vertices)
    bags = fill(UInt16(N+1), N)

    bag = 0
    for (i, vi) in enumerate(I)
        v = vertices[vi]
        d = degrees[v]
        queue[i] = Pair(v, d)
        keys[v] = i

        if d >= bag
            bags[bag+1:d+1] .= i
            bag = d + 1
        end
    end

    MMDQueue(queue, keys, bags, 0x0001)
end

"""
    pop!(q::MMDQueue)

Return the element with highest priority and increment all
bag boundaries to point to the next highest priority element.
"""
function pop!(q::MMDQueue)
    p = q.queue[q.min]
    q.min += 0x0001
    q.bags[0x0001:p.second+0x0001] .= q.min
    p
end

#=
decrement! is implemented by swapping v, in bag d, with the vertex 
at the boundary of bag d. The degree of v is decremented and the 
boundary of bag d is incremented so that v is no longer in bag d 
and now in bag d - 1.

Example:

Initial state: 
               ---to swap --
               |           |
               |           |
degree: [ ][ ][ ][ ][ ][ ][d][ ][ ]
queue : [ ][ ][p][ ][ ][ ][v][ ][ ]
bags  : [ ][b][ ][ ][ ][ ][ ][ ][ ]
keys  : [ ][ ][ ][ ][i][ ][ ][ ][ ]
ind   :     d  b     v     i


Final state: (D = d - 1, B = b + 1)

                  --- new boundary of bag d
                  |
                  |
degree: [ ][ ][D][|][ ][ ][ ][ ][ ]
queue : [ ][ ][v][|][ ][ ][p][ ][ ]
bags  : [ ][B][ ][|][ ][ ][ ][ ][ ]
keys  : [ ][ ][ ][|][b][ ][ ][ ][ ]
ind   :     d  b  |  v     i
=#

"""
    decrement!(q::MMDQueue, v::UInt16)

Decrease the degree of v by 1 and update the queue to maintain
the correct priority for vertices and bag structure.

See comments above docstring for implementation details.
"""
function decrement!(q::MMDQueue, v::UInt16)
    i = q.keys[v]           # index of v.
    if i >= q.min
        d = q.queue[i].second   # the degree of v.
        b = q.bags[d + 0x0001]  # index of bag d boundary.

        # Move the pair at the boundary to where v was.
        p = q.queue[b]
        q.queue[i] = p
        q.keys[p.first] = i

        # Place new pair for v at the boundary and
        # increment the boundary position.
        q.queue[b] = v => d - 0x0001
        q.keys[v] = b
        q.bags[d + 0x0001] += 0x0001
    end

    nothing
end

###
### Base.show implementation for MMDQueues.
###

function Base.show(io::IO, q::MMDQueue)
    if length(q.queue) <= 32
        show_verbose(io, q)
    else
        show_compact(io, q)
    end
end

"""
Print some of the fields of the given MMDQueue.
"""
function show_compact(io::IO, q::MMDQueue)
    ks = [p.first for p in q.queue]
    ds = [p.second for p in q.queue]
    println(io, "degree: ", Int.(ds))
    println(io, "vertex: ", Int.(ks))
    println(io, "bag   : ", Int.(q.bags))
    println(io, "key   : ", Int.(q.keys))
end

"""
Draw a verbose image of the given MMDQueue.
"""
function show_verbose(io::IO, q::MMDQueue)
    println()
    bar = "  |  "
    space = "     "

    # Draw the bag structure of the queue.
    display_array = -1*ones(Int, length(q.queue))
    for (i, b) in enumerate(q.bags)
        display_array[b] = i - 1
    end

    Bl = ""
    Bu = ""
    for i in display_array
        if i == -1
            Bl *= space
            Bu *= space
        else
            Bl *= bar
            Bu *= "  $i  "
        end
    end

    println(io, Bu)
    println(io, Bl)

    # Draw the queue.
    Q = ""
    D = ""
    for i = 1:length(q.queue)
        v, d = q.queue[i]
        Q *= "[" * lpad(v, 3) * "]"
        D *= "[" * lpad(d, 3) * "]"
    end

    println(io, Q)
    println(io, D)

    # Draw the arrow pointing to the next item of the queue.
    m = q.min
    min_arrow = space^(m - 1) * "  ^  "
    min_bar   = space^(m - 1) * bar
    min_base  = space^(m - 1) * " min "

    println(io, min_arrow)
    println(io, min_bar)
    println(io, min_base)
end

end