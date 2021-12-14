module Graphs

import LightGraphs as lg

export Graph, vertices, neighbours, neighbourhood, degree, lg
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

    adj_matrix::BitMatrix # true -> edge doesn't exists, false -> edge does exists
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

    # Initialise the graphs adjacency matrix.
    adj_matrix = trues(num_vertices, num_vertices)
    for v = 1:num_vertices
        for u in adj_list[v][1:degree[v]]
            adj_matrix[v, u] = false
            adj_matrix[u, v] = false
        end
    end

    # Create the graph.
    Graph(num_vertices, vertices, degree, adj_list, adj_matrix)
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
vertices(g::Graph) = @view g.vertices[1:g.num_vertices]
neighbours(g::Graph, v::Integer) = @view g.adj_list[v][1:g.degree[v]]

"""Return the adjacency matrix of the neighbourhood of v."""
function neighbourhood(g::Graph, v::Integer)
    N = neighbours(g, v)
    @view g.adj_matrix[N, N]
end

"""Remove a vertex from the graph."""
function remove_vertex!(g::Graph, v::UInt16)::Nothing
    # This assumes v is a vertex of g. who knows what happens if it isn't.

    # Find where v is in the vertices array.
    i = find_index(g.vertices, v)

    # Swap v with the n-th vertex and decrease n by 1.
    g.vertices[i] = g.vertices[g.num_vertices]
    g.vertices[g.num_vertices] = v
    g.num_vertices -= 1

    # Remove v from the neighbourhood of its neighbours.
    for n in neighbours(g, v)
        N = neighbours(g, n)
        i = find_index(N, v)

        N[i] = N[end]
        g.degree[n] -= 1

        # Remove edges the adjacency matrix
        # TODO: only use a triangular matrix.
        g.adj_matrix[v, n] = true
        g.adj_matrix[n, v] = true
    end

    nothing
end

"""Add an edge to the grapg G connecting u and v."""
function add_edge!(g::Graph, u::UInt16, v::UInt16)::Nothing
    # Increment the degrees of u and v.
    g.degree[u] += 1
    g.degree[v] += 1

    # Add v as a neighbour of u and vice versa.
    g.adj_list[u][g.degree[u]] = v
    g.adj_list[v][g.degree[v]] = u

    # Update adjacency matrix.
    g.adj_matrix[u, v] = false
    g.adj_matrix[v, u] = false

    nothing
end

"""Remove the dge connecting vertices u and v from the graph G."""
function remove_edge!(g::Graph, u::UInt16, v::UInt16)::Nothing
    # Update adjacency list for u.
    i = find_index(g.adj_list[u], v)
    g.adj_list[u][i] = g.adj_list[u][g.degree[u]]
    g.degree[u] -= 1

    # Update adjacency list for v.
    i = find_index(g.adj_list[v], u)
    g.adj_list[v][i] = g.adj_list[v][g.degree[v]]
    g.degree[v] -= 1

    # Update adjacency matrix.
    g.adj_matrix[u, v] = true
    g.adj_matrix[v, u] = true

    nothing
end

"""Get the number of edges that need to be added to make the neighbourhood a simplex."""
function simplexity(nbhd::SubArray)::Int64
    count = 0
    n = size(nbhd)[1]
    for i = 1:n-1
        for j = i+1:n
            @inbounds count += nbhd[i, j]
        end
    end
    count
end

"""Turn the neighbourhood of v into a simplex."""
function make_simplicial!(g::Graph, v::UInt16)::Vector{Tuple{UInt16, UInt16}}

    Nbhd = neighbourhood(g, v)
    N = neighbours(g, v)

    # TODO: use a static array here?
    edges_added = Vector{Tuple{UInt16, UInt16}}(undef, simplexity(Nbhd))

    d = g.degree[v]
    e = 1
    for i = 1:d-1
        for j = i+1:d
            if Nbhd[i, j]
                @inbounds ui = N[i]; uj = N[j]
                @inbounds edges_added[e] = (ui, uj)
                e += 1
                add_edge!(g, ui, uj)
            end
        end
    end

    edges_added
end

"""Remove the vertex v from the graph g and make it simplicial."""
function eliminate!(g::Graph, v::UInt16)::Vector{Tuple{UInt16, UInt16}}
    remove_vertex!(g, v)
    make_simplicial!(g, v)
end

"""
Add the last vertex to be eliminated back to g and restore its neighbourhood.

'edges_to_remove' is assumed to contain the edges added to g to make the vertex
simplicial when being eliminated.
"""
function restore_last_eliminated!(g::Graph, edges_to_remove::Vector{Tuple{UInt16, UInt16}})::Nothing
    g.num_vertices += 1
    v = g.vertices[g.num_vertices]

    for n in neighbours(g, v)
        g.degree[n] += 1
        g.adj_list[n][g.degree[n]] = v

        g.adj_matrix[v, n] = false
        g.adj_matrix[n, v] = false
    end

    for (ui, uj) in edges_to_remove
        remove_edge!(g, ui, uj)
    end
    
    nothing
end

###
### Utility functions
###

"""Get the index of a value v in an array A"""
function find_index(A::AbstractVector{T}, v::T)::Int64 where T
    i = 1
    while A[i] != v
        i += 1
    end
    i
end

end