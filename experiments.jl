using Statistics
include("activelearning.jl")
include("ida.jl")
include("genData.jl")

function runsingle(mtd, strat)
    # path to this directory -> CHANGE THIS
    dirpath = "/xyz/abc/counting-with-applications"
    # path to the installation of dct-policy (https://github.com/csquires/dct-policy) -> CHANGE THIS
    dctpath = "/xyz/abc/dct-policy"

    pathoutfile = "data/activelearning/" * mtd * "-" * strat * ".ans"
    outfile = open(pathoutfile, "a")
    println("Results are written to " * pathoutfile)
    gs = []
    if mtd == "hairball_plus_rand"
        gs = [100, 200, 300, 400, 500]
    end
    if mtd == "shanmugam_small_rand"
        gs = [10, 15, 20, 25, 30, 35, 40]
    end

    rep = 1
    
    for od in gs
        results = []
        times = []
        println(outfile, "Number of vertices: " * string(od) * " -----------")
        for nr in 1:rep
            infile = mtd * "-" * string(od) * "-" * string(rep) * ".gr"
            if !isfile("data/activelearning/" * mtd * "/dags/" * infile)
                continue
            end
            G = readgraph("data/activelearning/" * mtd * "/dags/" * infile)
            file = dirpath * "/data/activelearning/" * mtd * "/dags/" * infile
            print(outfile, file * ": ")
            println(file)
            # dct
            if strat == "dct"
                runscript = dctpath * "/run.sh"
                time = @elapsed res = parse(Int, split(read(`$runscript $file`, String), "\n")[1])
            end
            # coloring
            if strat == "col"
                runscript = dctpath * "/run_col.sh"
                time = @elapsed res = parse(Int, split(read(`$runscript $file`, String), "\n")[1])
            end
            # optsingle
            if strat == "os"
                time = @elapsed iset = launcher(G, optsingle)
                if !checker(G, iset)
                    println("ERROR: iset is not sufficient")
                    return
                end
                res = length(iset)
            end
            # entropy
            if strat == "ent"
                time = @elapsed iset = launcher(G, entropy)
                if !checker(G, iset)
                    println("ERROR: iset is not sufficient")
                    return
                end
                res = length(iset)
            end
            # minmax
            if strat == "minmax"
                time = @elapsed iset = launcher(G, minmax)
                if !checker(G, iset)
                    println("ERROR: iset is not sufficient")
                    return
                end
                res = length(iset)
            end

            if time > 1800 # stop if one case takes too long, could be improved
                close(outfile)
                return
            end
            
            println(outfile, "#int: " * string(res) * ", time in s: " * string(time))
            push!(times, time)
            push!(results, res)
            flush(outfile)
        end
        
        println(outfile, "Times: " * "n=" * string(od) * " mean=" * string(mean(times)) * " std=" * string(std(times)))
        println(outfile, "Results: " * "n=" * string(od) * " mean=" * string(mean(results)) * " std=" * string(std(results)))
        flush(outfile)
    end
    close(outfile)
end

function IDAexperiments()
    pathoutfile = "data/ida/results.ans"
    outfile = open(pathoutfile, "a")
    println("Results are written to " * pathoutfile)
    nn = [5, 10, 15, 30, 50, 100]
    for n in nn
        ltm = []
        mtm = []
        println(outfile, "Working on n = " * string(n) * ": ")
        rep = 1
        while length(mtm) < 100
            println(string(n) * " " * string(rep))
            dt = generateidadata(n, 3, 1000)
            df = DataFrame(dt, :auto)
            x = rand(1:n)
            y = rand(1:n)
            while y == x
                y = rand(1:n)
            end
            G = pcalg(df, 0.01, gausscitest)
            U = copy(G)
            U.ne = 0
            for i = 1:n
                filter!(j->has_edge(G, j, i), U.fadjlist[i])
                filter!(j->has_edge(G, i, j), U.badjlist[i])
                U.ne += length(U.fadjlist[i])
            end
            skip = false
            for component in connected_components(U)
                cc = induced_subgraph(G, component)
                for u in vertices(cc[1])
                    for v in inneighbors(cc[1], u)
                        if !has_edge(cc[1], u, v)
                            println("Not induced")
                            skip = true
                        end
                    end
                end
                if !ischordal(cc[1])
                    println("Undirected connected components are NOT chordal...Abort")
                    println("Are you sure the graph is a CPDAG?")
                    skip = true
                end
            end

            rep += 1

            # if the graph is not a CPDAG (we just check for chordality of undir. comp.),
            # we draw another graph
            if skip
                continue
            end

            # don't include init time
            # but use other variables
            # else glm stores result?!
            # because it is abnormally fast
            multiplicities(G, ((x+1) % n) + 1)
            loc_ida(G, ((x+1) % n) + 1, ((y+1) % n) + 1, df)
            
            tmp = @elapsed begin
                mults = multiplicities(G, x)
            end
            push!(mtm, tmp)

            tmp = @elapsed begin
                effects = loc_ida(G, x, y, df)
            end
            push!(ltm, tmp)
            
            if length(effects) != length(mults) # sanity check
                println("Something went wrong")
                println(effects)
                println(mults)
            end
        end
        
        println(outfile, "  [mult] avg time in s: " * string(mean(mtm)))
        println(outfile, "  [mult] std dev for time in s: " * string(std(mtm)))
        println(outfile, "  [local] avg time in s: " * string(mean(ltm)))
        println(outfile, "  [local] std dev for time in s: " * string(std(ltm)))
        flush(outfile)
    end

    df = DataFrame(CSV.File("data/ida/riboflavin.csv"))
    df = select(df, Not(:Column1))
    rename!(df, [Symbol("x$i") for i in 1:4089])
    G = pcalg(df, 0.01, gausscitest)

    # check that resulting graph has chordal induced components
    # the Clique-Picking functions rely on this fact
    n = nv(G)
    U = copy(G)
    U.ne = 0
    for i = 1:n
        filter!(j->has_edge(G, j, i), U.fadjlist[i])
        filter!(j->has_edge(G, i, j), U.badjlist[i])
        U.ne += length(U.fadjlist[i])
    end
    for component in connected_components(U)
        cc = induced_subgraph(G, component)
        for u in vertices(cc[1])
            for v in inneighbors(cc[1], u)
                if !has_edge(cc[1], u, v)
                    println("Not induced")
                end
            end
        end
        if !ischordal(cc[1])
            println("Undirected connected components are NOT chordal...Abort")
            println("Are you sure the graph is a CPDAG?")
        end
    end
    
    ltm = []
    mtm = []
    # don't include init time
    multiplicities(G, 2)
    loc_ida(G, 2, 3, df)
    
    for i = 2:4089
        println("mults for covariate " * string(i))
        tmp = @elapsed multiplicities(G, i)
        push!(mtm, tmp)
    end

    for i = 2:4089
        println("effects for covariate " * string(i))
        tmp = @elapsed loc_ida(G, i, 1, df)
        push!(ltm, tmp)
    end
    println(outfile, "riboflavin:")
    println(outfile, "  [mult] avg time in s: " * string(mean(mtm)))
    println(outfile, "  [mult] std dev for time in s: " * string(std(mtm)))
    println(outfile, "  [local] avg time in s: " * string(mean(ltm)))
    println(outfile, "  [local] std dev for time in s: " * string(std(ltm)))
    close(outfile)
end

function sampling_exp(mtd)
    pathoutfile = "data/sampling/" * mtd * ".ans"
    outfile = open(pathoutfile, "a")
    println("Results are written to " * pathoutfile)
    gs = [16, 32, 64, 128, 256, 512, 1024, 2048, 4096]
    rep = 1
    srep = 10
    for od in gs        
        cnttimes = []
        pretimes = []
        stimes = []
        println(outfile, "Number of vertices: " * string(od) * " -----------")

        # we make one initial run, which is not measured as Julia
        # compiles code during the first run
        for nr in 0:rep
            infile = mtd * "-n=" * string(od) * "-nr=" * string((nr < 1 || nr > rep) ? 1 : nr) * ".gr"
            if !isfile("data/sampling/" * mtd * "/chordalgraphs/" * infile)
                continue
            end
            G = readgraph("data/sampling/" * mtd * "/chordalgraphs/" * infile, true)
            file = "data/sampling/" * mtd * "/chordalgraphs/" * infile
            print(outfile, file * ": ")
            println(file)

            cnttime = 0.0
            pretime = 0.0
            
            if rand(1:2) == 1
                cnttime = @elapsed res = MECsize(G)
                # garbage collector is very lazy, has to be triggered
                # to properly clean up
                GC.gc()
                GC.gc()
                GC.gc()
                GC.gc()
                pretime = @elapsed pre = precomputation(G)
                GC.gc()
                GC.gc()
                GC.gc()
                GC.gc()
            else
                pretime = @elapsed pre = precomputation(G)
                GC.gc()
                GC.gc()
                GC.gc()
                GC.gc()
                cnttime = @elapsed res = MECsize(G)
                GC.gc()
                GC.gc()
                GC.gc()
                GC.gc()
            end

            ittime = 0
            
            for si in 1:srep
                ttime = @elapsed sampleDAG(G, pre)
                ittime += ttime
            end

            GC.gc()
            GC.gc()
            GC.gc()
            GC.gc()
            
            println(outfile, "cnttime: " * string(cnttime) * " pretime: " * string(pretime) * " average sampling time: " * string(ittime/srep))
            flush(outfile)
            
            if nr > 0 && nr <= rep
                push!(cnttimes, cnttime)
                push!(pretimes, pretime)
                push!(stimes, ittime/srep)
            end
        end

        println(outfile, "Counting: " * "n=" * string(od) * " mean=" * string(mean(cnttimes)) * "  std=" * string(std(cnttimes)))
        println(outfile, "Pre: " * "n=" * string(od) * " mean=" * string(mean(pretimes)) * "  std=" * string(std(pretimes)))
        println(outfile, "Sampling: " * "n=" * string(od) * " mean=" * string(mean(stimes)) * "  std=" * string(std(stimes)))
        flush(outfile)
    end
    close(outfile)
end

#####Active Learning experiments:
for mtd in ["shanmugam_small_rand", "hairball_plus_rand"]
    for strat in ["dct", "col", "os", "ent", "minmax"]
        runsingle(mtd, strat)
    end
end

##### IDA experiments:
IDAexperiments()

###### Sampling experiments:
for mtd in ["subtree-logn", "interval"]
     sampling_exp(mtd)
end
