using Random
using Graphs
using Distributions
using CausalInference
using DataFrames
using StatsModels
using GLM
using CSV
using DataStructures

include("activelearning.jl")

function randgraph(n, d, dist = "er")
    G = SimpleDiGraph(n)
    if dist == "er"
        p = d / (n-1)
        for a = 1:n, b = (a+1):n
            if rand() <= p
                add_edge!(G, a, b)
                add_edge!(G, b, a)
            end
        end
    end
    return G
end

function randDAG!(G)
    n = nv(G)
    ts = randperm(n)
    makeDAG!(G, ts)
end

function makeDAG!(G, ts)
    n = nv(G)
    for a in 1:n, b in inneighbors(G, a)
        if ts[b] < ts[a]
            rem_edge!(G, a, b)
        end
    end
end

function assignweights(G)
    n = nv(G)
    wg = [Dict{Integer, AbstractFloat}() for i = 1:n]
    ud = Uniform(1, 2)
    for a in vertices(G), b in outneighbors(G, a)
        wg[a][b] = rand(ud)
    end
    return wg
end

function sampledata(wg, s, ts)
    n = size(wg, 1)
    nd = Normal()
    dt = rand(nd, s, n)
    its = zeros(Int64, n)
    for i in 1:n
        its[ts[i]] = i
    end
    # assumes ts 1..n
    for a in its, b in wg[a]
        dt[1:s, b.first] += dt[1:s, a] * b.second
    end
    return dt
end

function generateidadata(n, d, s)
    G = randgraph(n, d, "er")
    ts = randperm(n)
    makeDAG!(G, ts)
    wg = assignweights(G)
    return sampledata(wg, s, ts)
end
