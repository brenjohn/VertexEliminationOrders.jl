export find_treewidth_from_order
export find_similar_groups

###
### functions to analyse elimination orders.
###

"""
    find_treewidth_from_order(G::AbstractGraph, order::Array{Symbol, 1})
Return the treewidth of `G` with respect to the elimination order in `order`.
"""
function find_treewidth_from_order(G::AbstractGraph, order::Array{Int, 1})
    G = deepcopy(G)
    labels = collect(1:nv(G))
    τ = 0
    for v_label in order
        v = findfirst(vl -> vl == v_label, labels)
        τ = max(τ, degree(G, v))
        eliminate!(G, labels, v)

        lb = mmd_plus(G)
        if lb >= 7
            println("Oops: lb is ", lb)
        end
    end
    τ
end

###
### Functions to analyse a graph
###

"""Return a vector containing the similar groups of G as defined by Yaun_2011"""
function find_similar_groups(G::SimpleGraph{Int})
    groups = Vector{Int}[]
    for v in vertices(G)
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

function are_similar(G::SimpleGraph, u::Int, v::Int)
    return Set(setdiff(all_neighbors(G, u), [v])) == Set(setdiff(all_neighbors(G, v), [u]))
end