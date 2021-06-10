using LyapunovExponents 
using Plots 
using LightGraphs 
using LinearAlgebra 

ϵ = 1000. 
E = -ϵ * laplacian_matrix(star_graph(5)) |> collect
P = [1 0 0; 0 0 0; 0 0 0 ]
net = Network(Lorenz, E, P)
eigenmodes(net)
