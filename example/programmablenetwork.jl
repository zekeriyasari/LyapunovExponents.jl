using LinearAlgebra: length
# Design of Programmable Networks 

using LyapunovExponents 
using LinearAlgebra
using LightGraphs
using GraphPlot
using DifferentialEquations
using Plots 
theme(:default) 
default(legend=:topleft, label=nothing) 
mt = (:circle, 5) 
ms = (:star5, 5) 
plotgraph(graph) = gplot(graph, nodelabel=1:nv(graph))
getwave(x, i, d) = getindex.(x, (i - 1) * d + 1) 

# Choose a topology 
n = 6 
graph = SimpleGraph(n) 
add_edge!(graph, 1, 2) 
add_edge!(graph, 1, 3) 
add_edge!(graph, 2, 3) 
add_edge!(graph, 4, 5) 
add_edge!(graph, 4, 6) 
add_edge!(graph, 5, 6) 
add_edge!(graph, 1, 4) 
plotgraph(graph) |> display
Ξ = collect(laplacian_matrix(graph)) .|> float

# Define the symmetry 
Q = [
    1 0 0 0 0 0; 
    0 0 1 0 0 0; 
    0 1 0 0 0 0; 
    0 0 0 1 0 0;
    0 0 0 0 0 1; 
    0 0 0 0 1 0
]
M = eigvecs(Q) 
G = inv(M) * Ξ * M 
B = G[1 : 2, 1 : 2] 
D = G[3 : 6, 3 : 6] 
λB = eigvals(B) 
λD = eigvals(D) 
ploteigs = plot() 
scatter!(λB, zeros(length(λB)), marker=mt, label="λt")
scatter!(λD, zeros(length(λD)), marker=ms, label="λs")
display(ploteigs)

# Note that this symmetry is stable. Choose another symmetry
n = 6 
graph = SimpleGraph(n) 
add_edge!(graph, 2, 5) 
add_edge!(graph, 1, 4) 
add_edge!(graph, 3, 6) 
add_edge!(graph, 1, 2) 
add_edge!(graph, 1, 3) 
add_edge!(graph, 4, 5) 
add_edge!(graph, 4, 6) 
plotgraph(graph) |> display
Ξ = collect(laplacian_matrix(graph)) .|> float
Q = [
    0 0 0 1 0 0; 
    0 0 0 0 1 0; 
    0 0 0 0 0 1; 
    1 0 0 0 0 0;
    0 1 0 0 0 0; 
    0 0 1 0 0 0
]
M = eigvecs(Q) 
G = inv(M) * Ξ * M 
B = G[1 : 2, 1 : 2] 
D = G[3 : 6, 3 : 6] 
λB = eigvals(B) 
λD = eigvals(D) 
ploteigs = plot() 
scatter!(λB, zeros(length(λB)), marker=mt, label="λt")
scatter!(λD, zeros(length(λD)), marker=ms, label="λs")
display(ploteigs)

# This mode is also stable 