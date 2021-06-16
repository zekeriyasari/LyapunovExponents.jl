# Programming network

using LyapunovExponents 
using LinearAlgebra
using LightGraphs
using GraphPlot
using Plots 
theme(:default) 
default(legend=:topleft, label=nothing) 
mt = (:circle, 5) 
ms = (:star5, 5) 
plotgraph(graph) = gplot(graph, nodelabel=1:nv(graph), layout=circular_layout)
getwave(x, i, d) = getindex.(x, (i - 1) * d + 1)  

function scale(W0)
    W = copy(W0)
    for i in 1 : size(W, 1)
        W[i, :] ./=  W[i, i]
    end
    W
end

function sweepweight(w, W, M, i, j)
    λB = Vector{Vector{Float64}}(undef, length(w))
    λD = Vector{Vector{Float64}}(undef, length(w))
    for (l, wi) in enumerate(w)
        # W[2, 5] = W[5, 2] =  -wi
        # W[2, 2] -= sum(W[2, :])
        # W[5, 5] -= sum(W[5, :])
        W[i, j] = W[j, i] = -wi
        W[i, i] = -sum(W[i, [k for k ∈ 1 : size(W, 2) if k ≠ i]])
        W[j, j] = -sum(W[j, [k for k ∈ 1 : size(W, 2) if k ≠ j]])
        Ξ = scale(W)
        G = inv(M) * Ξ * M
        B = G[1 : 3, 1 : 3]
        D = G[4 : 6, 4 : 6]
        λB[l] = eigvals(B)
        λD[l] = eigvals(D)
    end
    λB, λD
end

# Define a graph and check eigenvalue distribution 
n = 6 
graph = cycle_graph(n)
add_edge!(graph, 1, 4)
add_edge!(graph, 2, 6)
add_edge!(graph, 3, 5)
plotgraph(graph) |> display
W = collect(laplacian_matrix(graph)) .|> float
Ξ = scale(W)

# Choose a symmetry 
Q = [
    0 0 0 1 0 0; 
    0 0 1 0 0 0;
    0 1 0 0 0 0;
    1 0 0 0 0 0;
    0 0 0 0 0 1;
    0 0 0 0 1 0;
]   
M = eigvecs(Q) 
G = inv(M) * Ξ * M 
B = G[1:3, 1:3] 
D = G[4:6, 4:6]
λB = eigvals(B) 
λD = eigvals(D) 
ploteigs = plot() 
scatter!(λB, zeros(length(λB)), marker=mt)
scatter!(λD, zeros(length(λD)), marker=ms)

# Sweep the graph 
w = range(0, 4, length=41)
λB, λD = sweepweight(w, W, M, 1, 4)
ploteig = plot()
plot!(w, getindex.(λB, 1), marker=mt,  label="λ1t")
plot!(w, getindex.(λB, 2), marker=mt,  label="λ2t")
plot!(w, getindex.(λB, 3), marker=mt,  label="λ3t")
plot!(w, getindex.(λD, 1), marker=ms,  label="λ1s")
plot!(w, getindex.(λD, 2), marker=ms,  label="λ2s")
plot!(w, getindex.(λD, 3), marker=ms,  label="λ3s")
xlabel!("w")

# Change the link 
wi = 3.
i, j = 1, 4 
W[i, j] = W[j, i] = -wi
W[i, i] = -sum(W[i, [k for k ∈ 1 : size(W, 2) if k ≠ i]])
W[j, j] = -sum(W[j, [k for k ∈ 1 : size(W, 2) if k ≠ j]])
Ξ = scale(W)
G = inv(M) * Ξ * M
B = G[1 : 3, 1 : 3]
D = G[4 : 6, 4 : 6]
λB = eigvals(B)
λD = eigvals(D)
ploteigs = plot() 
scatter!(λB, zeros(length(λB)), marker=mt)
scatter!(λD, zeros(length(λD)), marker=ms)
