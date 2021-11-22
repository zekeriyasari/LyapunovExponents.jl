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

# To determine the stability of the synchronization, we need to find the threshold point for MSF 

ds = Lorenz(σ=10, β=8/3, ρ=35) 
P = diagm([1., 0., 0.])
Ψ(η) = msf(ds, η, P)
ηr = range(0, 20, length=21) 
Λr = Ψ.(ηr) 
ηtr = ηr[findfirst(Λr .< 0)]    # Threshold for η
plotmsf = plot() 
plot!(ηr, Λr, marker=(:circle, 3), xlabel="η", label="Λ")
hline!([0.], linestyle=:dash, linewidth=2) 
display(plotmsf) 


let 
    # Construct a graph 
    n = 8 
    graph = cycle_graph(n) 
    add_edge!(graph, 2, 7)
    add_edge!(graph, 3, 6)
    W = collect(laplacian_matrix(graph)) .|> float 
    Ξ = scale(copy(W))
    plotgraph(graph) |> display

    # ----------------------------- Q1 symmetry -------------- #

    # Define a network symmetry 
    Q = [                   # Symmetry S1
        0 0 0 0 0 0 0 1;
        0 0 0 0 0 0 1 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 1 0 0 0;
        0 0 0 1 0 0 0 0;
        0 0 1 0 0 0 0 0;
        0 1 0 0 0 0 0 0;
        1 0 0 0 0 0 0 0;
    ]
    M = eigvecs(Q)

    # Q symmetry is not stable. Let us sweep the link 2 and 5 
    i, j = 1, 8 
    ws = range(0, 10, length=51) 
    λsB = Vector{Vector{Float64}}(undef, length(ws))
    λsD = Vector{Vector{Float64}}(undef, length(ws))
    for (l, wi) in enumerate(ws)
        W[i, j] = W[j, i] = -wi 
        W[i, i] = -sum(W[i, [k for k ∈ 1 : size(W, 2) if k ≠ i]])
        W[j, j] = -sum(W[j, [k for k ∈ 1 : size(W, 2) if k ≠ j]])
        Ξ = scale(W) 
        G = inv(M) * Ξ * M
        B = G[1:4, 1:4]
        D = G[5:8, 5:8]
        λsB[l] = eigvals(B) 
        λsD[l] = eigvals(D) 
    end
    ploteigsweep = plot(xlabel="w$i$j", title="S1 Symmetry") 
    plot!(ws, getindex.(λsB, 1), marker=mt, label="λt1")
    # plot!(ws, getindex.(λsB, 2), marker=mt, label="λt2")
    # plot!(ws, getindex.(λsD, 1), marker=ms, label="λs1")
    plot!(ws, getindex.(λsD, 2), marker=ms, label="λs2")
    # plot!(ws, getindex.(λsD, 3), marker=ms, label="λs3")
    # plot!(ws, getindex.(λsD, 4), marker=ms, label="λs4")
    display(ploteigsweep)

    # Define a network symmetry
    wt = 1
    W[i, j] = W[j, i] = -wt 
    W[i, i] = -sum(W[i, [k for k ∈ 1 : size(W, 2) if k ≠ i]])
    W[j, j] = -sum(W[j, [k for k ∈ 1 : size(W, 2) if k ≠ j]])
    Ξ = scale(W) 
    G = inv(M) * Ξ * M
    B = G[1:4, 1:4]
    D = G[5:8, 5:8]
    λB = eigvals(B) 
    λD = eigvals(D) 
    ploteig = plot() 
    scatter!(λB, zeros(length(λB)), marker=mt, label="λt")
    scatter!(λD, zeros(length(λD)), marker=ms, label="λs")
    display(ploteig)
    λt = λB[1] 
    λs = λD[2] 
    ϵr1 = ηtr ./ [λt, λs]    # Couplin strength range 
    @show ϵr1
    ϵ1 = sum(ϵr1) / 2

    # ----------------------------- Q2 symmetry -------------- #

    # Define a network symmetry 
    W = collect(laplacian_matrix(graph)) .|> float 
    Ξ = scale(copy(W))
    Q = [                   # Symmetry S1
        0 0 0 1 0 0 0 0;
        0 0 1 0 0 0 0 0;
        0 1 0 0 0 0 0 0;
        1 0 0 0 0 0 0 0;
        0 0 0 0 0 0 0 1;
        0 0 0 0 0 0 1 0;
        0 0 0 0 0 1 0 0;
        0 0 0 0 1 0 0 0;
    ]
    M = eigvecs(Q)

    # Q symmetry is not stable. Let us sweep the link 2 and 5 
    i, j = 5, 8 
    ws = range(0, 10, length=51) 
    λsB = Vector{Vector{Float64}}(undef, length(ws))
    λsD = Vector{Vector{Float64}}(undef, length(ws))
    for (l, wi) in enumerate(ws)
        W[i, j] = W[j, i] = -wi 
        W[i, i] = -sum(W[i, [k for k ∈ 1 : size(W, 2) if k ≠ i]])
        W[j, j] = -sum(W[j, [k for k ∈ 1 : size(W, 2) if k ≠ j]])
        Ξ = scale(W) 
        G = inv(M) * Ξ * M
        B = G[1:4, 1:4]
        D = G[5:8, 5:8]
        λsB[l] = eigvals(B) 
        λsD[l] = eigvals(D) 
    end
    ploteigsweep = plot(xlabel="w$i$j", title="S2 Symmetry") 
    plot!(ws, getindex.(λsB, 1), marker=mt, label="λt1")
    # plot!(ws, getindex.(λsB, 2), marker=mt, label="λt2")
    # plot!(ws, getindex.(λsB, 3), marker=mt, label="λt3")
    # plot!(ws, getindex.(λsD, 1), marker=ms, label="λs1")
    plot!(ws, getindex.(λsD, 2), marker=ms, label="λs2")
    # plot!(ws, getindex.(λsD, 3), marker=ms, label="λs3")
    display(ploteigsweep)

    # Define a network symmetry
    wt = 5
    W[i, j] = W[j, i] = -wt 
    W[i, i] = -sum(W[i, [k for k ∈ 1 : size(W, 2) if k ≠ i]])
    W[j, j] = -sum(W[j, [k for k ∈ 1 : size(W, 2) if k ≠ j]])
    Ξ = scale(W) 
    G = inv(M) * Ξ * M
    B = G[1:4, 1:4]
    D = G[5:8, 5:8]
    λB = eigvals(B) 
    λD = eigvals(D) 
    ploteig = plot() 
    scatter!(λB, zeros(length(λB)), marker=mt, label="λt")
    scatter!(λD, zeros(length(λD)), marker=ms, label="λs")
    display(ploteig)
    λt = λB[1] 
    λs = λD[2] 
    ϵr2 = ηtr ./ [λt, λs]    # Couplin strength range 
    @show ϵr2
    ϵ2 = sum(ϵr2) / 2
end 
