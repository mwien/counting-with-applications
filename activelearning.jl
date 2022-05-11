using Graphs
include("counting.jl")
include("sampling.jl")

function find_parents(G, mapping, trueDAG, v)
    parents = Vector{Int}()
    for x in inneighbors(G, v)
        if has_edge(trueDAG, mapping[x], mapping[v]) 
            push!(parents, x)
        end
    end
    return parents
end

function remaining_components(G, v, parents)
    n = nv(G)
    # do a graph search from v without edges to parents
    visited = falses(n)
    for x in parents
        visited[x] = true
    end
    D = Vector{Int}()
    push!(D, v)
    visited[v] = true
    i = 1
    while(i <= length(D))
        x = D[i]
        for y in inneighbors(G, x)
            if !visited[y]
                visited[y] = true
                push!(D, y)
            end
        end
        i += 1
    end
    # call subproblems function on P \cup D with clique P \cup v
    subgraph, mp = induced_subgraph(G, vcat(parents, D))
    imp = zeros(Int64, n)
    for i in 1:length(mp)
        imp[mp[i]] = i
    end
    K = Set{Int}(map(x -> imp[x], parents))
    push!(K, imp[v])
    sps = subproblems(subgraph, K)
    sps = map(x -> mp[x], sps)
    
    # add A \cup P as subproblem and return
    !isempty(parents) && push!(sps, setdiff(vertices(G), D))
    return sps
end

function count_intsize(cc, v, parents, memo, fmemo)
    G = cc[1] # graph
    mapping = cc[2] # mapping to original vertex numbers
    n = nv(G)
    size = 1
    for H in remaining_components(G, v, parents)
        HH = induced_subgraph(G, H)
        size *= count((HH[1], map(x -> mapping[x], HH[2])), memo, fmemo)
    end
    return size
end

function cliques(G, nb)
    cl = Vector{Vector{Int64}}() 
    push!(cl, Vector{Int64}())
    # do this more elegantly
    for v in nb
        ncl = Vector{Vector{Int64}}()
        for x in cl
            isc = true
            for y in x
                if !has_edge(G, v, y)
                    isc = false
                    break
                end
            end
            if isc
                nx = copy(x)
                push!(nx, v)
                push!(ncl, nx)
            end
        end
        cl = cat(cl, ncl, dims = 1)
    end
    return cl
end

function minmax(cc, memo, fmemo)
    G = cc[1] # graph
    mapping = cc[2] # mapping to original vertex numbers
    n = nv(G)
    mn = -1
    mni = -1
    for v in vertices(G)
        mx = 0
        for parents in cliques(G, inneighbors(G, v))
            mx = max(mx, count_intsize(cc, v, parents, memo, fmemo))         
        end
        if mn == -1 || mx < mn
            mn = mx
            mni = v
        end
    end
    return mni
end

function entropy(cc, memo, fmemo)
    G = cc[1] # graph
    mapping = cc[2] # mapping to original vertex numbers
    n = nv(G)
    mx = 0
    mxi = -1
    for v in vertices(G)
        ent = 0
        L = count(cc, memo, fmemo)
        for parents in cliques(G, inneighbors(G, v))
            x = count_intsize(cc, v, parents, memo, fmemo)
            ent -= x / L * log2(x / L)
        end
        if mx < abs(ent)
            mx = abs(ent)
            mxi = v
        end
    end
    return mxi
end

function optsingle(cc, memo, fmemo) 
    G = cc[1] # graph
    mapping = cc[2] # mapping to original vertex numbers
    n = nv(G)
    mn = -1
    mni = -1
    for v in vertices(G)
        mx = 0
        for parents in cliques(G, inneighbors(G, v))
            sumedges = 0
            for H in remaining_components(G, v, parents)
                HH = induced_subgraph(G, H)
                sumedges += ne(HH[1])
            end
            mx = max(mx, sumedges)
        end
        if mn == -1 || mx < mn
            mn = mx
            mni = v
        end
    end
    return mni
end

function actlearn_simulator(cc, inttarget::Function, trueDAG, memo, fmemo)
    G = cc[1] # graph
    mapping = cc[2] # mapping to original vertex numbers
    if nv(G) == 1
        return Set()
    end
    v = inttarget(cc, memo, fmemo)
    interventions = Set()
    push!(interventions, mapping[v])
    parents = find_parents(G, mapping, trueDAG, v)
    for H in remaining_components(G, v, parents)
        HH = induced_subgraph(G, H)
        union!(interventions, actlearn_simulator((HH[1], map(x -> mapping[x], HH[2])), inttarget, trueDAG, memo, fmemo))
    end
    return interventions
end

function launcher(trueDAG, inttarget) # make this more general
    n = nv(trueDAG)
    G = SimpleDiGraph(n)
    for v in vertices(trueDAG)
        for x in outneighbors(trueDAG, v) 
            add_edge!(G, v, x)
            add_edge!(G, x, v)
        end 
    end
    memo = Dict{Set, BigInt}() #mapping set of vertices -> AMO sum
    fmemo = zeros(BigInt, n)
    U = copy(G)
    U.ne = 0
    for i = 1:n
        filter!(j->has_edge(G, j, i), U.fadjlist[i])
        filter!(j->has_edge(G, i, j), U.badjlist[i])
        U.ne += length(U.fadjlist[i])
    end
    tinterventions = Set()
    for component in connected_components(U)
        cc = induced_subgraph(U, component)
        if !ischordal(cc[1])
            println("Undirected connected components are NOT chordal...Abort")
            println("Are you sure the graph is a CPDAG?")
            # is there anything more clever than just returning?
            return
        end
        union!(tinterventions, actlearn_simulator(cc, inttarget, trueDAG, memo, fmemo))
    end

    return tinterventions
end

function meek!(G)
    n = nv(G)
    done = false
    while !done
        done = true
        # rule 1:
        # a -> b - c => a -> b -> c
        for a = 1:n, b in outneighbors(G, a)
            if !has_edge(G, b, a)
                for c in outneighbors(G, b)
                    if has_edge(G, c, b) && !has_edge(G, a, c) && !has_edge(G, c, a)
                        rem_edge!(G, c, b)
                        done = false
                    end
                end
            end
        end
        # rule 2:
        # a -> b -> c and a - c => a -> c
        for a = 1:n, b in outneighbors(G, a)
            if !has_edge(G, b, a)
                for c in outneighbors(G, b)
                    if !has_edge(G, c, b) && has_edge(G, a, c) && has_edge(G, c, a)
                        rem_edge!(G, c, a)
                        done = false
                    end
                end
            end
        end
        # rule 3 is not necessary
    end
end

function checker(trueDAG, iset)
    n = nv(trueDAG)
    G = SimpleDiGraph(n)
    for v in vertices(trueDAG)
        if v in iset
            for x in outneighbors(trueDAG, v)
                add_edge!(G, v, x)
            end
        else
            for x in outneighbors(trueDAG, v)
                add_edge!(G, v, x)
                if !(x in iset)
                    add_edge!(G, x, v)
                end
            end
        end
    end
    meek!(G)
    for v in vertices(G)
        for x in outneighbors(G, v)
            if has_edge(G, x, v)
                return false
            end
        end
    end
    return true
end
