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

function runsim1(ϵ)
    n = 8 
    P = diagm([1., 0., 0.])
    graph = cycle_graph(n) 
    add_edge!(graph, 2, 7)
    add_edge!(graph, 3, 6)
    plotgraph(graph) |> display
    W = collect(laplacian_matrix(graph)) .|> float 
    Ξ = scale(copy(W))
    nodes = [Lorenz(σ=10, β=8/3, ρ=35) for i in 1 : n]
    net = Network(nodes, -ϵ * Ξ, P)
    sol = solvenet(net, (0., 1000.))
    t, x = sol.t, sol.u
    d = dimension(net.nodes[1])
    default(label=nothing)
    plotres = plot(layout=(n,n), size=(1500, 600)) 
    idx = 1
    for i in 1 : n 
        for j in 1 : n
            plot!(t, abs.(getwave(x, i, d) - getwave(x, j, d)), subplot=idx)
            idx += 1
        end
    end
    plotres
end 

function runsim2(ϵ)
    n = 8 
    P = diagm([1., 0., 0.])
    graph = cycle_graph(n) 
    add_edge!(graph, 2, 7)
    add_edge!(graph, 3, 6)
    add_edge!(graph, 5, 8)
    plotgraph(graph) |> display
    W = collect(laplacian_matrix(graph)) .|> float 
    wt = 5
    i, j = 5, 8
    W[i, j] = W[j, i] = -wt 
    W[i, i] = -sum(W[i, [k for k ∈ 1 : size(W, 2) if k ≠ i]])
    W[j, j] = -sum(W[j, [k for k ∈ 1 : size(W, 2) if k ≠ j]])
    Ξ = scale(copy(W))
    nodes = [Lorenz(σ=10, β=8/3, ρ=35) for i in 1 : n]
    net = Network(nodes, -ϵ * Ξ, P)
    sol = solvenet(net, (0., 1000.))
    t, x = sol.t, sol.u
    d = dimension(net.nodes[1])
    default(label=nothing)
    plotres = plot(layout=(n,n), size=(1500, 600)) 
    idx = 1
    for i in 1 : n 
        for j in 1 : n
            plot!(t, abs.(getwave(x, i, d) - getwave(x, j, d)), subplot=idx)
            idx += 1
        end
    end
    plotres
end 

runsim1(15) 
runsim2(15)

