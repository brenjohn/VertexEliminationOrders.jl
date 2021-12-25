module Graphs

#=
This module defines a custom graph struct, tailored for a depth first 
search through the set of all vertex elimination orders of a graph.
It aims to provide graph operations which are both central to such a 
search algorithm and are non-allocating.

The implemented interface allows one to view the current vertices of a 
graph, the neighbours of a vertex and and its degree. It also enables 
one to eliminate a vertex and restore the previous vertex which was 
eliminated.
=#

#=
Info on the implementation:

The Graph struct contains the following fields:

    num_vertices - num vertices remaining in the graph.
    vertices     - a vector to hold remaining and removed vertices.
    degree       - a vector to hold the degrees of vertices.
    adj_list     - an adjacency list
    adj_memory   - a matrix to track which and when edges are added.
    num_updates  - num updates applied to the graph to get to the 
                   current state. Initialisation is the first update,
                   the rest are eliminations.

In a depth first search of vertex elimination orders, after 
initialisation, various vertices are eliminated from the graph and 
restored, resulting in various intermediate states.

At any given state, the vertices that remain in the graph are stored
in the first 'num_vertices' elements of the vector 'vertices'. The
other elements of this vector maintain the order in which the 
eliminated vertices were removed.

The connectivity of the remaining vertices is described by the
adjacency list. Namely, the neighbours of a vertex v in the graph are
given by 'adj_list[v][1:degree[v]]'. If v is a vertex which was 
eliminated, then this expression gives the neighbours v had when it
was eliminated.

When eliminating a vertex, two things happen to facilitate the 
restoration of previously eliminated vertices:

    1) a counter is incremented to track the number of updates made 
       to the graph to get to the current state.

    2) elements of a matrix, called adjacency memory matrix, are updated
       to record which edges were added to make the vertex simplicial and
       which update these edges were added in.

The adjacency memory matrix is similar to an adjacency matrix but not
equivalent. It not only contains information regarding the current connectivity
of the graph but also temporal information about when certain edges were added
to the graph. If 'M = adj_memory' and both 'i' and 'j' are remaining vertices
in the graph, then 'M[i, j] == 0' means there is no edge between 'i' and 'j' 
while 'M[i, j] == u', where 'u > 0', means there is an edge between them and
it was aded during update 'u'. Note, the initialisation of the graph is taken
as update number 1. If either 'i' or 'j' are eliminated, then 'M[i,j]' 
indicates the state of the associated edge when the earliest of these updates
happened.

Safety checks, such making sure a vertex isnt already removed before eliminating,
have been ommitted in the interest of performance and such assumptions are 
assumed to be guaranteed by the depth first search implementation.
=#

import LightGraphs as lg

export Graph, num_vertices, vertices, neighbours, degree, lg
export eliminate!, restore_last_eliminated!

###
### Graph struct
###

"""struct to represent a graph"""
mutable struct Graph
    num_vertices::UInt16
    vertices::Vector{UInt16}

    degree::Vector{UInt16}
    adj_list::Vector{Vector{UInt16}}

    adj_memory::Matrix{UInt16}
    num_updates::UInt16
end

"""Graph Constructor"""
function Graph(G::lg.SimpleGraph)
    num_vertices = lg.nv(G)
    vertices = lg.vertices(G)
    degree = lg.degree(G)

    # Initialise the graphs adjacency list.
    adj_list = [Array{UInt16, 1}(undef, num_vertices-1) for i = 1:num_vertices]
    for v = 1:num_vertices
        adj_list[v][1:degree[v]] = lg.all_neighbors(G, v)
    end

    # Initialise the graphs adjacency memory.
    adj_memory = zeros(UInt16, num_vertices, num_vertices)
    for v = 1:num_vertices
        for u in adj_list[v][1:degree[v]]
            adj_memory[v, u] = 0x0001
            adj_memory[u, v] = 0x0001
        end
    end

    # Create the graph.
    Graph(num_vertices, vertices, degree, adj_list, adj_memory, 0x0001)
end

function Base.show(io::IO, g::Graph)
    compact = get(io, :compact, false)
    if !compact
        println(io, "The graph has ", g.num_vertices, " nodes")
        println(io, "and ", Int(sum(g.degree[g.vertices[1:g.num_vertices]])/2), " edges.")
    end
end

###
### Graph interface functions.
###

degree(g::Graph, v::Integer) = g.degree[v]
num_vertices(g::Graph) = g.num_vertices
vertices(g::Graph) = @view g.vertices[1:g.num_vertices]
neighbours(g::Graph, v::Integer) = @view g.adj_list[v][1:g.degree[v]]

"""Return the adjacency matrix of the neighbourhood of v."""
function neighbourhood(g::Graph, v::Integer)
    N = neighbours(g, v)
    @view g.adj_memory[N, N]
end

"""Remove a vertex from the graph."""
function remove_vertex!(g::Graph, v::UInt16)::Nothing
    # This assumes v is a vertex of g. who knows what happens if it isn't.

    # Find where v is in the vertices array.
    i = find_index(g.vertices, v)

    # Swap v with the n-th vertex and decrease n by 1.
    g.vertices[i] = g.vertices[g.num_vertices]
    g.vertices[g.num_vertices] = v
    g.num_vertices -= 0x0001

    # Remove v from the neighbourhood of its neighbours.
    for n in neighbours(g, v)
        N = neighbours(g, n)
        i = find_index(N, v)

        N[i] = N[end]
        g.degree[n] -= 0x0001
    end

    nothing
end

"""Add an edge to the grapg G connecting u and v."""
function add_edge!(g::Graph, u::UInt16, v::UInt16)::Nothing
    # Increment the degrees of u and v.
    g.degree[u] += 0x0001
    g.degree[v] += 0x0001

    # Add v as a neighbour of u and vice versa.
    g.adj_list[u][g.degree[u]] = v
    g.adj_list[v][g.degree[v]] = u

    # Update adjacency matrix.
    g.adj_memory[u, v] = g.num_updates
    g.adj_memory[v, u] = g.num_updates

    nothing
end

"""Remove the dge connecting vertices u and v from the graph G."""
function remove_edge!(g::Graph, u::UInt16, v::UInt16)::Nothing
    # Update adjacency list for u.
    i = find_index(g.adj_list[u], v)
    g.adj_list[u][i] = g.adj_list[u][g.degree[u]]
    g.degree[u] -= 0x0001

    # Update adjacency list for v.
    i = find_index(g.adj_list[v], u)
    g.adj_list[v][i] = g.adj_list[v][g.degree[v]]
    g.degree[v] -= 0x0001

    # Update adjacency matrix.
    g.adj_memory[u, v] = 0x0000
    g.adj_memory[v, u] = 0x0000

    nothing
end

"""Get the number of edges that need to be added to make the neighbourhood a simplex."""
function simplexity(nbhd::SubArray)::Int64
    count = 0x0000
    n = size(nbhd)[1]
    for i = 1:n-1
        for j = i+1:n
            @inbounds if nbhd[i, j] == 0x0000
                count += 0x0001
            end
        end
    end
    count
end

"""
Turn the neighbourhood of v into a simplex.
The number of edges added to do so is also returned.
"""
function make_simplicial!(g::Graph, v::UInt16)::UInt16
    Nbhd = neighbourhood(g, v)
    N = neighbours(g, v)

    d = g.degree[v]
    e = 0
    for i = 1:d-1
        for j = i+1:d
            if Nbhd[i, j] == 0x0000
                @inbounds ui = N[i]; uj = N[j]
                add_edge!(g, ui, uj)
                e += 0x0001
            end
        end
    end

    e
end

"""Remove the vertex v from the graph g and make it simplicial."""
function eliminate!(g::Graph, v::UInt16)::UInt16
    g.num_updates += 0x0001
    remove_vertex!(g, v)
    make_simplicial!(g, v)
end

"""
Add the last vertex to be eliminated back to g and restore its neighbourhood.

'edges_to_remove' is assumed to contain the edges added to g to make the vertex
simplicial when being eliminated.
"""
function restore_last_eliminated!(g::Graph)::Nothing
    g.num_vertices += 0x0001
    v = g.vertices[g.num_vertices]
    N = neighbours(g, v)

    for n in N
        g.degree[n] += 0x0001
        g.adj_list[n][g.degree[n]] = v
    end

    Nbhd = neighbourhood(g, v)
    d = g.degree[v]
    for i = 1:d-1
        for j = i+1:d
            if Nbhd[i, j] == g.num_updates
                @inbounds ui = N[i]; uj = N[j]
                remove_edge!(g, ui, uj)
            end
        end
    end
    
    g.num_updates -= 0x0001

    nothing
end

###
### Utility functions
###

"""Get the index of a value v in an array A"""
function find_index(A::AbstractVector{T}, v::T)::Int64 where T
    i = 1
    @inbounds while A[i] != v
        i += 1
    end
    i
end

end