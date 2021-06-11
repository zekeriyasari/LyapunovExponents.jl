# This file includes stability modes of a network 
using LyapunovExponents
using LightGraphs 
using GraphPlot
using Plots 
using LinearAlgebra

# Construct graph 
gr = cycle_graph(5) 
add_edge!(gr, 1, 3) 
add_edge!(gr, 1, 4) 
add_edge!(gr, 2, 5) 
gplot(gr, layout=circular_layout, nodelabel=1:nv(gr))

# Construct dynamic network 
ϵ = 50
E = collect(laplacian_matrix(gr))
P = [1 0 0; 0 0 0; 0 0 0 ]
net = Network(Lorenz, -ϵ * E, P)

nsteps = Int(3e4)
ntrsteps = Int(1e4) 
dt = 0.001 
λ = range(0, 100, length=51) 
Λ = map(λi -> λi => msf(net.nodes[1], λi, net.P), λ)
@info "Done calculation"
plot(λ, getfield.(Λ, :second), marker=(:circle, 3))
hline!([0], linestyle=:dash)


# # Solve the network 
# sol = solvenet(net, (0., 50.))
# t, x = sol.t, sol.u

# # Plot solutions 
# plot(layout = (3, 1))
# plot!(t, abs.(getindex.(x, 4) - getindex.(x, 13)), subplot=1)
# plot!(t, abs.(getindex.(x, 7) - getindex.(x, 10)), subplot=2)
# plot!(t, abs.(getindex.(x, 1) - getindex.(x, 4)), subplot=3)


# Define permutatin matrix 
Q = [ 
    1 0 0 0 0; 
    0 0 0 0 1; 
    0 0 0 1 0; 
    0 0 1 0 0; 
    0 1 0 0 0
]
M = eigvecs(Q)
G = inv(M) * E * M
B = G[1 : 2, 1 : 2] 
D = G[3 : 5, 3 : 5] 
valB =  map(((k, λi),) -> λi => msf(net.nodes[k], λi, net.P), enumerate(eigvals(B)))
valD =  map(((k, λi),) -> λi => msf(net.nodes[k], λi, net.P), enumerate(eigvals(D)))
@show valB 
@show valD 
nothing
