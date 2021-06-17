# # Topological Control of Networks 

# Load packages and define some defaults 
using LyapunovExponents
using LightGraphs
using GraphPlot
using LinearAlgebra
using Plots 
theme(:default) 
default(legend = :topleft)
mt = (:circle, 5) 
ms = (:star5, 5) 

getwave(x, i, d) = getindex.(x, (i - 1) * d + 1)

plotgraph(graph::AbstractGraph) = gplot(graph, nodelabel=1:nv(graph), layout=circular_layout)
plotgraph(net::Network) = plotgraph(SimpleDiGraph(net.E))

function scale(mat) 
    W = copy(mat) 
    for i in 1 : size(mat, 2)
        W[i, :] ./= W[i, i]
    end
    W
end 

let 
    # Construct a graph 
    n = 6 
    graph = cycle_graph(n) 
    add_edge!(graph, 1, 4)
    add_edge!(graph, 2, 6)
    add_edge!(graph, 3, 5)
    W = collect(laplacian_matrix(graph)) .|> float 
    Ξ = scale(copy(W))
    plotgraph(graph) |> display

    # ----------------------------- Q1 symmetry -------------- #

    # Define a network symmetry 
    Q = [                   # Symmetry S1
        1 0 0 0 0 0;
        0 0 0 0 0 1;
        0 0 0 0 1 0;
        0 0 0 1 0 0;
        0 0 1 0 0 0;
        0 1 0 0 0 0;
    ]
    M = eigvecs(Q)

    # Q symmetry is not stable. Let us sweep the link 2 and 5 
    i, j = 2, 3
    ws = range(0, 5, length=51) 
    λsB = Vector{Vector{Float64}}(undef, length(ws))
    λsD = Vector{Vector{Float64}}(undef, length(ws))
    for (l, wi) in enumerate(ws)
        W[i, j] = W[j, i] = -wi 
        W[i, i] = -sum(W[i, [k for k ∈ 1 : size(W, 2) if k ≠ i]])
        W[j, j] = -sum(W[j, [k for k ∈ 1 : size(W, 2) if k ≠ j]])
        Ξ = scale(W) 
        G = inv(M) * Ξ * M
        B = G[1:2, 1:2]
        D = G[3:6, 3:6]
        λsB[l] = eigvals(B) 
        λsD[l] = eigvals(D) 
    end
    ploteigsweep = plot(xlabel="w", title="S1 Symmetry") 
    plot!(ws, getindex.(λsB, 1), marker=mt, label="λt1")
    # plot!(ws, getindex.(λsB, 2), marker=mt, label="λt2")
    # plot!(ws, getindex.(λsD, 1), marker=ms, label="λs1")
    plot!(ws, getindex.(λsD, 2), marker=ms, label="λs2")
    # plot!(ws, getindex.(λsD, 3), marker=ms, label="λs3")
    # plot!(ws, getindex.(λsD, 4), marker=ms, label="λs4")
    display(ploteigsweep)

    # ----------------------------- Q2 symmetry -------------- #

    # Define a network symmetry 
    Q = [               # Symmetry S2
        0 0 0 1 0 0;
        0 0 1 0 0 0;
        0 1 0 0 0 0;
        1 0 0 0 0 0;
        0 0 0 0 0 1;
        0 0 0 0 1 0;
    ]
    M = eigvecs(Q)

    # Q symmetry is not stable. Let us sweep the link 2 and 5 
    i, j = 2, 3
    ws = range(0, 5, length=51) 
    λsB = Vector{Vector{Float64}}(undef, length(ws))
    λsD = Vector{Vector{Float64}}(undef, length(ws))
    for (l, wi) in enumerate(ws)
        W[i, j] = W[j, i] = -wi 
        W[i, i] = -sum(W[i, [k for k ∈ 1 : size(W, 2) if k ≠ i]])
        W[j, j] = -sum(W[j, [k for k ∈ 1 : size(W, 2) if k ≠ j]])
        Ξ = scale(W) 
        G = inv(M) * Ξ * M
        B = G[1:3, 1:3]
        D = G[4:6, 4:6]
        λsB[l] = eigvals(B) 
        λsD[l] = eigvals(D) 
    end
    ploteigsweep = plot(xlabel="w", title="S2 Symmetry") 
    plot!(ws, getindex.(λsB, 1), marker=mt, label="λt1")
    # plot!(ws, getindex.(λsB, 2), marker=mt, label="λt2")
    # plot!(ws, getindex.(λsB, 3), marker=mt, label="λt3")
    # plot!(ws, getindex.(λsD, 1), marker=ms, label="λs1")
    plot!(ws, getindex.(λsD, 2), marker=ms, label="λs2")
    # plot!(ws, getindex.(λsD, 3), marker=ms, label="λs3")
    display(ploteigsweep)
end 
