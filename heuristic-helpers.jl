function label(g)
   #returns digraph w/ labels
    gOut = MetaDiGraph(nv(g))
    for i = 1:nv(g) #copy node labels over
       set_props!(gOut,i,props(g,i)) 
    end
    unlab = Dict() #dict counting # unlabelled neighbors for node i
    weights = Dict()
    queue = []
    if nv(g) == 1
        set_props!(gOut,1,Dict(:weight=>1))
        return gOut
    end
    for i = 1:nv(g)
        unlab[i] = neighbors(g,i)
        weights[i] = 0
        if length(neighbors(g,i)) == 1
            append!(queue,i)
        end
    end
    while length(queue) > 0
        curr = popfirst!(queue)
        nw = [weights[i] for i in neighbors(g,curr)]
        if length(unlab[curr]) == 1
            weights[curr] = sum(nw) + 1
            add_edge!(gOut,convert(Int64,unlab[curr][1]),curr)
        else
            weights[curr] = sum(nw) + 1 - findmax(nw)[1]
        end
        for n in convert(Array{Int64,1},neighbors(g,curr))
            unlab[n] = filter(x->x.!=curr,unlab[n])
            if length(unlab[n]) == 1
                push!(queue,n)
            end
        end
        set_props!(gOut,curr,Dict(:weight=>weights[curr]))
    end
    root = findfirst([length(inneighbors(gOut,i)) for i in 1:nv(gOut)].==0)
    queue = [root]
    set_props!(gOut,root,Dict(:layer=>1))
    while length(queue) > 0
        v = popfirst!(queue)
        for child in outneighbors(gOut,v)
            set_prop!(gOut,child,:layer,props(gOut,v)[:layer]+1)
            push!(queue,child)
        end
    end
    return gOut
end

function plotmg(g,prop=:name)
    if prop == :none
        layout=(args...)->spring_layout(args...; C=5)
        gplot(g,layout=layout)
    else
        gplot(g,nodelabel=[props(g,i)[prop] for i in 1:nv(g)])
    end
end

function numAdded(assem,parent)
    added = 0
    if !(parent in assem.singles) && !(parent in assem.threads)
        added += 1
        for child in outneighbors(assem.graph,parent)
            added += numAdded(assem,child)
        end
    end
    return added
end

function fast_heuristic_run(s,R,O,tstop)
    states = [s]
    T = 1
    while !isempty(s) && T<tstop
        sp = heuristic_policy(s,R,O)
        [nv(sp[i].graph) for i in 1:length(sp)]
        push!(states,sp)
        s = sp
        T+=1
    end
    return (states,T)
end

function greedyPolicy(const_state,R,S)
   initAssign(const_state[1],0,R)
    sp = partitionAssem(const_state[1])
    return sp
end

function greedy_heuristic_run(s,R,O,tstop)
    states = [s]
    T = 1
    while !isempty(s) && T<tstop
        sp = greedyPolicy(s,R,O)
        push!(states,sp)
        s = sp
        T+=1
    end
    return (states,T)
end