include("activelearning.jl")

function multiplicities(G, x)
    n = nv(G)
    memo = Dict{Set, BigInt}()
    fmemo = zeros(BigInt, n)
    mults = Dict{Vector{Int64}, BigInt}()

    bp = Vector{Int64}()
    for p in inneighbors(G, x)
        if !has_edge(G, x, p)
            push!(bp, p)
        end
    end
    
    U = copy(G)
    U.ne = 0
    for i = 1:n
        filter!(j->has_edge(G, j, i), U.fadjlist[i])
        filter!(j->has_edge(G, i, j), U.badjlist[i])
        U.ne += length(U.fadjlist[i])
    end
    base = 1
    xc = undef
    for component in connected_components(U)
        if !(x in component)
            cc = induced_subgraph(U, component)
            # if !ischordal(cc[1])
            #     println("Undirected connected components are NOT chordal...Abort")
            #     println("Are you sure the graph is a CPDAG?")
            #     # is there anything more clever than just returning?
            #     return
            # end
            base *= count(cc, memo, fmemo)
        else
            xc = component
        end
    end
    cc = induced_subgraph(U, xc)
    H = cc[1]
    mp = cc[2]
    ix = -1
    for i in 1:length(mp)
        if mp[i] == x
            ix = i
            break
        end
    end
    for parents in cliques(H, inneighbors(H, ix))
        parentset = Vector{Int64}()
        for z in parents
            push!(parentset, mp[z])
        end
        for z in bp
            push!(parentset, z)
        end
        mults[parentset] = base * count_intsize(cc, ix, parents, memo, fmemo)
    end
    return mults
end

function loc_ida(G, x, y, df)
    n = nv(G)
    effects = Dict{Vector{Int64}, Float64}()
    bp = Vector{Int64}()
    up = Vector{Int64}()
    for p in inneighbors(G, x)
        if !has_edge(G, x, p)
            push!(bp, p)
        else
            push!(up, p)
        end
    end

    for parents in cliques(G, up)
        resp = Term(Meta.parse("x" * string(y)))
        pr = ["x" * string(x)]
        parentset = Vector{Int64}()
        for z in parents
            push!(pr, "x" * string(z))
            push!(parentset, z)
        end
        for z in bp
            push!(pr, "x" * string(z))
            push!(parentset, z)
        end
        pr = map(x -> Meta.parse(x), pr)
        pred = Tuple(Term.(pr))
        effects[parentset] = coef(lm(FormulaTerm(resp, pred), df))[2]
    end
    return effects
end
