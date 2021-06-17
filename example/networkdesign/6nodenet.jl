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

# # To determine the stability of the synchronization, we need to find the threshold point for MSF 

# ds = Lorenz(σ=10, β=8/3, ρ=35) 
# P = [
#     0 0 0;
#     1 0 0;
#     0 0 0
#     ]
# Ψ(η) = msf(ds, η, P)
# ηr = range(0, 20, length=41) 
# Λr = Ψ.(ηr) 
# idx = findfirst(Λr .< 0)
# ηtr = sum(ηr[idx - 1 : idx]) / 2    # Threshold for η
# plotmsf = plot() 
# plot!(ηr, Λr, marker=(:circle, 3), xlabel="η", label="Λ")
# hline!([0.], linestyle=:dash, linewidth=2) 
# display(plotmsf) 

# Construct a graph 
n = 6 
graph = cycle_graph(n) 
add_edge!(graph, 1, 4)
add_edge!(graph, 2, 6)
add_edge!(graph, 3, 5)
W = collect(laplacian_matrix(graph)) .|> float 
wt = 0.5
i, j = 1, 4
W[i, j] = W[j, i] = -wt 
W[i, i] = -sum(W[i, [k for k ∈ 1 : size(W, 2) if k ≠ i]])
W[j, j] = -sum(W[j, [k for k ∈ 1 : size(W, 2) if k ≠ j]])
Ξ = scale(copy(W))
plotgraph(graph) |> display

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
G = inv(M) * Ξ * M
B = G[1:2, 1:2]
D = G[3:6, 3:6]
λB = eigvals(B) 
λD = eigvals(D) 
ploteig = plot() 
scatter!(λB, zeros(length(λB)), marker=mt, label="λt")
scatter!(λD, zeros(length(λD)), marker=ms, label="λs")
display(ploteig)
λt = λB[1] 
λs = λD[2] 
ϵr1 = ηtr ./ [λt, λs]    # Coupling strength range 
@show ϵr1
ϵ1 = sum(ϵr1) / 2

# Simulate the networks 
ϵ = 10.5
ϵr1[1] ≤ ϵ ≤ ϵr1[2] || error("The coupling is not in the range")
nodes = [Lorenz(σ=10, β=8/3, ρ=35) for i in 1 : n]
net = Network(nodes, -ϵ * Ξ, P)
sol = solvenet(net, (0., 1000.))
t, x = sol.t, sol.u
d = dimension(net.nodes[1])
plotres = plot(layout=(4,2))
plot!(t, abs.(getwave(x, 2, d) - getwave(x, 6, d)), label="e26", subplot=1)
plot!(t, abs.(getwave(x, 3, d) - getwave(x, 5, d)), label="e35", subplot=2)
plot!(t, abs.(getwave(x, 2, d) - getwave(x, 3, d)), label="e23", subplot=3)
plot!(t, abs.(getwave(x, 5, d) - getwave(x, 6, d)), label="e56", subplot=4)
plot!(t, abs.(getwave(x, 1, d) - getwave(x, 4, d)), label="e14", subplot=5)
plot!(t, abs.(getwave(x, 1, d) - getwave(x, 2, d)), label="e12", subplot=6)
plot!(t, abs.(getwave(x, 4, d) - getwave(x, 5, d)), label="e45", subplot=7)
display(plotres)

# --------------------------------------------------------------------------

# Define a network symmetry
wt = 5
i, j = 1, 4
W[i, j] = W[j, i] = -wt 
W[i, i] = -sum(W[i, [k for k ∈ 1 : size(W, 2) if k ≠ i]])
W[j, j] = -sum(W[j, [k for k ∈ 1 : size(W, 2) if k ≠ j]])
Ξ = scale(copy(W))
Q = [               # Symmetry S2
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
ploteig = plot() 
scatter!(λB, zeros(length(λB)), marker=mt, label="λt")
scatter!(λD, zeros(length(λD)), marker=ms, label="λs")
display(ploteig)
λt = λB[1] 
λs = λD[2] 
ϵr2 = ηtr ./ [λt, λs]    # Couplin strength range 
@show ϵr2
ϵ2 = sum(ϵr2) / 2

# Simulate the networks 
ϵ = 10.5
ϵr2[1] ≤ ϵ ≤ ϵr2[2] || error("The coupling is not in the range")
nodes = [Lorenz(σ=10, β=8/3, ρ=35) for i in 1 : n]
net = Network(nodes, -ϵ * Ξ, P)
sol = solvenet(net, (0., 1000.))
t, x = sol.t, sol.u
d = dimension(net.nodes[1])
plotres = plot(layout=(4,2))
plot!(t, abs.(getwave(x, 1, d) - getwave(x, 4, d)), label="e14", subplot=1)
plot!(t, abs.(getwave(x, 2, d) - getwave(x, 3, d)), label="e23", subplot=2)
plot!(t, abs.(getwave(x, 5, d) - getwave(x, 6, d)), label="e56", subplot=3)
plot!(t, abs.(getwave(x, 1, d) - getwave(x, 2, d)), label="e12", subplot=4)
plot!(t, abs.(getwave(x, 1, d) - getwave(x, 6, d)), label="e16", subplot=5)
plot!(t, abs.(getwave(x, 2, d) - getwave(x, 6, d)), label="e26", subplot=6)
plot!(t, abs.(getwave(x, 3, d) - getwave(x, 5, d)), label="e35", subplot=7)
display(plotres)



