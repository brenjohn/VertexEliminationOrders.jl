export Graph, neighbourhood, degree

mutable struct Graph
    num_vertices::UInt16
    vertices::Vector{UInt16}

    degree::Vector{UInt16}
    neighbours::Vector{Vector{UInt16}}

    simplexity::Vector{Vector{Tuple{UInt16, UInt16}}} # edges that need to be added to turn neighbourhood of vertex into a simplex.
end

function Graph(G::lg.SimpleGraph)
    num_vertices = lg.nv(G)
    vertices = lg.vertices(G)

    degree = lg.degree(G)

    neighbours = [Array{UInt16, 1}(undef, num_vertices-1) for i = 1:num_vertices]
    for v = 1:num_vertices
        neighbours[v][1:degree[v]] = lg.all_neighbors(G, v)
    end

    simplexity = [Tuple{UInt16, UInt16}[] for v = 1:num_vertices]
    for v = 1:num_vertices
        for ni = 1:degree[v]-1
            for nj = ni+1:degree[v]
                ui = neighbours[v][ni]
                uj = neighbours[v][nj]
                if !lg.has_edge(G, ui, uj)
                    edge = ui < uj ? (ui, uj) : (uj, ui)
                    push!(simplexity[v], edge)
                end
            end
        end
    end

    Graph(num_vertices, vertices, degree, neighbours, simplexity)
end

function Base.show(io::IO, g::Graph)
    compact = get(io, :compact, false)
    if !compact
        println(io, "The graph has ", g.num_vertices, " nodes")
        println(io, "and ", Int(sum(g.degree[g.vertices[1:g.num_vertices]])/2), " edges.")
    end
end

vertices(g::Graph) = @view g.vertices[1:g.num_vertices]
neighbourhood(g::Graph, v::Integer) = @view g.neighbours[v][1:g.degree[v]]
degree(g::Graph, v::Integer) = g.degree[v]

function remove_vertex!(g::Graph, v::UInt16)
    # This assumes v is a vertex of g. who knows what happens if it isn't.

    # Find where v is in the vertices array.
    i = 1
    while g.vertices[i] != v
        i += 1
    end

    # Swap v with the n-th vertex and decrease n by 1.
    g.vertices[i] = g.vertices[g.num_vertices]
    g.vertices[g.num_vertices] = v
    g.num_vertices -= 1

    # Remove v from the neighbourhood of its neighbours.
    for n in neighbourhood(g, v)
        i = 1
        N = neighbourhood(g, n)
        while N[i] != v
            i += 1
        end

        N[i] = N[end]
        g.degree[n] -= 1

        # Remove edges from simplexity of n if it contains v
        filter!(edge -> !(v in edge), g.simplexity[n])
    end

    nothing
end

function add_edge!(g::Graph, u::UInt16, v::UInt16)
    # This assumes u is less than v and that the edge (u, v) doesn't already exist.

    # Update the simplexity of common neighbours of u and v.
    Nv = neighbourhood(g, v)
    Nu = neighbourhood(g, u)
    for n in Nu
        if n in Nv
            println(u, " and ", v, " affecting ", n)
            # remove (u, v) from simplexity[n]
            i = 1
            while g.simplexity[n][i] != (u, v)
                i += 1
            end
            deleteat!(g.simplexity[n], i)
        end
    end

    for n in Nu
        if !(v in g.neighbours[n])
            edge = n < v ? (n, v) : (v, n)
            push!(g.simplexity[u], edge)
        end
    end

    for n in Nv
        if !(u in g.neighbours[n])
            edge = n < u ? (n, u) : (u, n)
            push!(g.simplexity[v], edge)
        end
    end

    # Increment the degrees of u and v.
    g.degree[u] += 1
    g.degree[v] += 1

    # Add v as a neighbour of u and vice versa.
    g.neighbours[u][g.degree[u]] = v
    g.neighbours[v][g.degree[v]] = u

    nothing
end

function remove_edge!(g::Graph, u::UInt16, v::UInt16)
    # This assumes u < v.
    i = 1
    while g.neighbours[u][i] != v
        i += 1
    end
    g.neighbours[u][i] = g.neighbours[u][g.degree[u]]
    g.degree[u] -= 1

    i = 1
    while g.neighbours[v][i] != u
        i += 1
    end
    g.neighbours[v][i] = g.neighbours[v][g.degree[v]]
    g.degree[v] -= 1

    # Update simplexity of neighbours.
    Nv = neighbourhood(g, v)
    Nu = neighbourhood(g, u)
    for n in Nu
        if n in Nv
            push!(g.simplexity[n], (u, v))
        end
    end

    for n in Nu
        filter!(e -> !(v in e), g.simplexity[n])
    end
    for n in Nv
        filter!(e -> !(u in e), g.simplexity[n])
    end

    nothing
end

function make_simplicial!(g::Graph, v::UInt16)
    for (ui, uj) in g.simplexity[v]
        println("adding ", ui, " and ", uj)
        add_edge!(g, ui, uj)
    end
    nothing
end

function eliminate!(g::Graph, v::UInt16)
    remove_vertex!(g, v)
    make_simplicial!(g, v)

    nothing
end

function restore_last_eliminated!(g::Graph)
    # I think the simplexity is incorrectly updated here
    # Should be an update for every edge connecting v
    g.num_vertices += 1
    v = g.vertices[g.num_vertices]

    for n in neighbourhood(g, v)
        g.degree[n] += 1
        g.neighbours[n][g.degree[n]] = v
    end

    for (ui, uj) in g.simplexity[v]
        remove_edge!(g, ui, uj)
    end
end