export graph_from_gr, graph_to_gr

"""
    graph_from_gr(filename::String)

Read a graph from the provided gr file.
"""
function graph_from_gr(filename::String)
    lines = readlines(filename)

    # Create a Graph with the correct number of vertices.
    num_vertices = parse.(Int, split(lines[1], ' ')[3])
    G = SimpleGraph(num_vertices)

    # Add an edge to the graph for every other line in the file.
    for line in lines[2:end]
        src, dst = parse.(Int, split(line, ' '))
        add_edge!(G, src, dst)
    end

    G
end

"""
    graph_to_gr(G::AbstractGraph, filename::String)

Write the provided graph to a file in gr format.
"""
function graph_to_gr(G::AbstractGraph, filename::String)
    open(filename, "w") do io
        write(io, "p tw $(nv(G)) $(ne(G))\n")
        for e in edges(G)
            write(io, "$(e.src) $(e.dst)\n")
        end
    end
end