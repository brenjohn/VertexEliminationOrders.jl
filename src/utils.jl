export find_treewidth_from_order
export find_similar_groups
export graph_from_gr, graph_to_gr, load_example_graph


###
### IO functions for reading and writing graphs in gr format.
###

"""
    graph_from_gr(filename::String)

Read a graph from the provided gr file.
"""
function graph_from_gr(filename::String)
    lines = readlines(filename)

    # Create a Graph with the correct number of vertices.
    num_vertices = parse.(Int, split(lines[1], ' ')[3])
    G = lg.SimpleGraph(num_vertices)

    # Add an edge to the graph for every other line in the file.
    for line in lines[2:end]
        src, dst = parse.(Int, split(line, ' '))
        lg.add_edge!(G, src, dst)
    end

    G
end

"""
    graph_to_gr(G::AbstractGraph, filename::String)

Write the provided graph to a file in gr format.
"""
function graph_to_gr(G::lg.AbstractGraph, filename::String)
    open(filename, "w") do io
        write(io, "p tw $(nv(G)) $(ne(G))\n")
        for e in lg.edges(G)
            write(io, "$(e.src) $(e.dst)\n")
        end
    end
end

"""
    load_example_graph(graph_name::String="test")

Return one of the example graphs from the `VertexEliminationOrders` package.

# Possible arguments:
- bristlecone\\_48\\_1-16-1\\_0
- rectangular\\_2x2\\_1-2-1\\_0
- rectangular\\_4x4\\_1-16-1\\_0
- sycamore\\_53\\_8\\_0
- sycamore\\_53\\_20\\_0
- test
"""
function load_example_graph(graph_name::String="test")
    filename = joinpath(dirname(dirname(@__FILE__)), "examples", graph_name * ".gr")
    graph_from_gr(filename)
end


###
### functions to analyse elimination orders.
###

"""
    find_treewidth_from_order(G::AbstractGraph, order::Array{Symbol, 1})

Return the treewidth of `G` with respect to the elimination order in `order`.
"""
function find_treewidth_from_order(G::lg.AbstractGraph, order::Array{<:Integer, 1})
    G = deepcopy(G)
    labels = collect(1:lg.nv(G))
    τ = 0
    for v_label in order
        v = findfirst(vl -> vl == v_label, labels)
        τ = max(τ, lg.degree(G, v))
        lg_eliminate!(G, labels, v)
    end
    τ
end


###
### Functions to analyse a graph.
###

"""Return a vector containing the similar groups of G as defined by Yaun_2011"""
function find_similar_groups(G::lg.SimpleGraph{Int})
    groups = Vector{Int}[]
    for v in lg.vertices(G)
        v_in_a_group = false
        for group in groups
            if are_similar(G, group[1], v)
                push!(group, v)
                v_in_a_group = true
                break
            end
        end
        v_in_a_group || push!(groups, [v])
    end

    groups
end

function are_similar(G::lg.SimpleGraph, u::Int, v::Int)
    return Set(setdiff(lg.all_neighbors(G, u), [v])) == Set(setdiff(lg.all_neighbors(G, v), [u]))
end

"""
    cliqueness(G::AbstractGraph, v::Integer)
Return the number of edges that need to be added to `G` in order to make the neighborhood of 
vertex `v` a clique.
"""
function cliqueness(g::lg.AbstractGraph, v::Integer)::Int
    neighborhood = lg.all_neighbors(g, v)::Array{Int64, 1}
    count = 0
    for i in 1:length(neighborhood)-1
        for j in i+1:length(neighborhood)
            vi = neighborhood[i]
            ui = neighborhood[j]
            if !lg.has_edge(g, ui, vi)::Bool
                count += 1
            end
        end
    end
    count
end


###
### Functions to modify a graph.
###

"""
    eliminate!(graph, labels, vertex)

Connects the neighbours of the given vertex into a clique before removing 
it from the graph.

The array of vertex labels are also updated to reflect the 
reording of vertex indices when the graph is updated.
"""
function lg_eliminate!(g::lg.AbstractGraph, labels, v)
    Nᵥ = lg.all_neighbors(g, v)::Array{Int64, 1}
    for i = 1:length(Nᵥ)-1
        vi = Nᵥ[i]
        for j = i+1:length(Nᵥ)
            vj = Nᵥ[j]
            lg.add_edge!(g, vi, vj)
        end
    end

    lg.rem_vertex!(g, v)
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
function lg_eliminate!(g::lg.AbstractGraph, labels, c_map, v)
    Nᵥ = lg.all_neighbors(g, v)::Array{Int64, 1}
    for i = 1:length(Nᵥ)-1
        vi = Nᵥ[i]
        for j = i+1:length(Nᵥ)
            vj = Nᵥ[j]

            # Try add an edge connecting vi and ui. If successful, update `c_map`.
            edge_added = lg.add_edge!(g, vi, vj)
            if edge_added
                # Common neighbours of vi and vj have one less edge to add
                # when being eliminated after vi and ui are connected.
                Nvi = lg.all_neighbors(g, vi)::Array{Int64, 1}
                Nvj = lg.all_neighbors(g, vj)::Array{Int64, 1}
                for n in Nvi
                    if n in Nvj
                        c_map[n] -= 1
                    end
                end

                # ui and vi are now neighbours so their cliqueness may increase.
                for n in Nvi
                    if !(n == vj) && !(lg.has_edge(g, n, vj)::Bool)
                        c_map[vi] += 1
                    end
                end
                for n in Nvj
                    if !(n == vi) && !(lg.has_edge(g, n, vi)::Bool)
                        c_map[vj] += 1
                    end
                end
            end
        end
    end

    # Removing v from G means it's also removed from its neighbour's neighbourhood, so their 
    # cliqueness may be reduced.
    for n in Nᵥ
        Nₙ = lg.all_neighbors(g, n)::Array{Int64, 1}
        for u in Nₙ
            if !(u == v)
                if !lg.has_edge(g, v, u)
                    c_map[n] -= 1
                end
            end
        end
    end

    lg.rem_vertex!(g, v)
    c_map[v] = c_map[end]
    labels[v] = labels[end]
    g
end