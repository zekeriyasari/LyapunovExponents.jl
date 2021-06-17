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

# Construct a graph 
n = 5 
graph = cycle_graph(n) 
add_edge!(graph, 1, 3)
add_edge!(graph, 1, 4)
add_edge!(graph, 2, 5)
W = collect(laplacian_matrix(graph)) .|> float 
Ξ = scale(W)
plotgraph(graph) |> display

# Define a network symmetry
Q = [
    1 0 0 0 0;
    0 0 0 0 1; 
    0 0 0 1 0; 
    0 0 1 0 0; 
    0 1 0 0 0
]
M = eigvecs(Q)
G = inv(M) * Ξ * M
B = G[1:2, 1:2]
D = G[3:5, 3:5]
λB = eigvals(B) 
λD = eigvals(D) 
ploteig1 = plot() 
scatter!(λB, zeros(length(λB)), marker=mt)
scatter!(λD, zeros(length(λD)), marker=ms)

# Q symmetry is not stable. Let us sweep the link 2 and 5 
i, j = 2, 5
ws = range(0, 2, length=21) 
λsB = Vector{Vector{Float64}}(undef, length(ws))
λsD = Vector{Vector{Float64}}(undef, length(ws))
for (l, wi) in enumerate(ws)
    W[i, j] = W[j, i] = -wi 
    W[i, i] = -sum(W[i, [k for k ∈ 1 : size(W, 2) if k ≠ i]])
    W[j, j] = -sum(W[j, [k for k ∈ 1 : size(W, 2) if k ≠ j]])
    Ξ .= scale(W) 
    G .= inv(M) * Ξ * M
    B .= G[1:2, 1:2]
    D .= G[3:5, 3:5]
    λsB[l] = eigvals(B) 
    λsD[l] = eigvals(D) 
end
ploteigsweep = plot(xlabel="w") 
plot!(ws, getindex.(λsB, 1), marker=mt, label="λt")
plot!(ws, getindex.(λsB, 2), marker=mt, label="λt2")
plot!(ws, getindex.(λsD, 1), marker=ms, label="λs1")
plot!(ws, getindex.(λsD, 2), marker=ms, label="λs")
plot!(ws, getindex.(λsD, 3), marker=ms, label="λs3")
display(ploteigsweep)

# Note that when the weight between the nodes 2 and 5 exceeds 1, the symmetry becomes stable. Let us choose 1.5
wt = 1.5
i, j = 2, 5
W[i, j] = W[j, i] = -wt 
W[i, i] = -sum(W[i, [k for k ∈ 1 : size(W, 2) if k ≠ i]])
W[j, j] = -sum(W[j, [k for k ∈ 1 : size(W, 2) if k ≠ j]])
Ξ = scale(W)
G .= inv(M) * Ξ * M 
B .= G[1 : 2, 1 : 2] 
D .= G[3 : 5, 3 : 5] 
λB .= eigvals(B) 
λD .= eigvals(D) 
λt = λB[1] 
λs = λD[2] 
ϵr = ηtr ./ [λt, λs]    # Couplin strength range 
@show ϵr
ϵ1 = 9.6

# Construct a network 
nodes = [Lorenz(σ=10, β=8/3, ρ=35) for i in 1 : n]
net = Network(nodes, -ϵ1 * Ξ, P)
sol = solvenet(net, (0., 1000.))
t, x = sol.t, sol.u
d = dimension(net.nodes[1])
plotres = plot(layout=(4,1))
plot!(t, abs.(getwave(x, 2, d) - getwave(x, 5, d)), label="e25", subplot=1)
plot!(t, abs.(getwave(x, 3, d) - getwave(x, 4, d)), label="e34", subplot=2)
plot!(t, abs.(getwave(x, 1, d) - getwave(x, 2, d)), label="e12", subplot=3)
plot!(t, abs.(getwave(x, 1, d) - getwave(x, 3, d)), label="e13", subplot=4)
display(plotres)

