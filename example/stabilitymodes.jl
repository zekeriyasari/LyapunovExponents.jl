using LyapunovExponents 
using Plots 
using LightGraphs 
using GraphPlot
using LinearAlgebra 

# Consruct a graph
gr = SimpleGraph(5)
add_edge!(gr, 1, 2)

# Construt a network 
ϵ = 10. 
E = -ϵ * laplacian_matrix(star_graph(5)) |> collect
P = [1 0 0; 0 0 0; 0 0 0]
net = Network(Lorenz, E, P)

λ = eigvals(-net.E) 
Λ = map(((k, λi),) -> λi => msf(net.nodes[k], λi, P), enumerate(λ))

# Solve network 
sol = solvenet(net, (0., 100.), saveat=0.01)
t, x = sol.t, sol.u

plot(getindex.(x, 1), getindex.(x, 4))

