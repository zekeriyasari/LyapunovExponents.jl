# # Topological Control of Dynamical Networks 

# In this notebook, we will consider topological control of dynamical networks. 

# First load required packages and  define some helper functions 

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

# We can now choose a network graph 

n1 = 5 
graph1 = cycle_graph(n1) 
add_edge!(graph1, 1, 3) 
add_edge!(graph1, 1, 4) 
add_edge!(graph1, 2, 5) 
plotgraph(graph1) |> display
Ξ1 = collect(laplacian_matrix(graph1)) .|> float

# The cluster pattern can be associated with the symmetries of the network. Each symmetry can be associated with a
# permutation matrix. 

Q1 = [
    1 0 0 0 0; 
    0 0 0 0 1; 
    0 0 0 1 0; 
    0 0 1 0 0;
    0 1 0 0 0
]
M1 = eigvecs(Q1) 
G1 = inv(M1) * Ξ1 * M1 
B1 = G1[1 : 2, 1 : 2] 
D1 = G1[3 : 5, 3 : 5] 
ploteigs1 = plot() 
λB1 = eigvals(B1) 
λD1 = eigvals(D1) 
scatter!(λB1, zeros(length(λB1)), marker=mt, label="tr", xlabel="λ")
scatter!(λD1, zeros(length(λD1)), marker=ms, label="syn", xlabel="λ")
display(ploteigs1)

# Note that the eigeen value condtion is not satisfied since the first eigenvalue of the transversal subspace and least
# nontrivial eigenvalue of the synchronization subspace overlaps. To make them seperated, we can control a link of the
# network. Let the control link be the one between the node 2 and 5. 

wr1 = range(0, 2, length=21)
i, j = 2, 5 
λrB1 = Vector{Vector{Float64}}(undef, length(wr1))
λrD1 = Vector{Vector{Float64}}(undef, length(wr1))
for (k, wri) in enumerate(wr1) 
    Ξ1[i, j] = Ξ1[j, i] = -wri
    Ξ1[i, i] = -sum(Ξ1[i, [k for k ∈ 1 : size(Ξ1, 2) if k ≠ i]])
    Ξ1[j, j] = -sum(Ξ1[j, [k for k ∈ 1 : size(Ξ1, 2) if k ≠ j]])
    G1 .= inv(M1) * Ξ1 * M1 
    B1 .= G1[1 : 2, 1 : 2] 
    D1 .= G1[3 : 5, 3 : 5] 
    λrB1[k] = eigvals(B1) 
    λrD1[k] = eigvals(D1) 
end 
ploteigrange1 = plot(xlabel="w")
plot!(wr1,  getindex.(λrB1, 1), marker=mt, label="λt1")
plot!(wr1,  getindex.(λrB1, 2), marker=mt, label="λt2")
plot!(wr1,  getindex.(λrD1, 1), marker=ms, label="λs1")
plot!(wr1,  getindex.(λrD1, 2), marker=ms, label="λs2")
plot!(wr1,  getindex.(λrD1, 3), marker=ms, label="λs3")
display(ploteigrange1)

# Note that for $w > 1$, the eigenvalue corresponding to the least transversal mode becomes larger than the eigenvalue
# corresponding to the least nontrivial eigenmode of the synchronization subspace. Let us choose $ w = 2 $. And,
# determine the range of $\epsilon$, 

wt1 = 2 
Ξ1[i, j] = Ξ1[j, i] = -wt1
Ξ1[i, i] = -sum(Ξ1[i, [k for k ∈ 1 : size(Ξ1, 2) if k ≠ i]])
Ξ1[j, j] = -sum(Ξ1[j, [k for k ∈ 1 : size(Ξ1, 2) if k ≠ j]])
G1 .= inv(M1) * Ξ1 * M1 
B1 .= G1[1 : 2, 1 : 2] 
D1 .= G1[3 : 5, 3 : 5] 
ϵr1  = ηtr ./  [eigvals(B1)[1], eigvals(D1)[2]] 
ϵ1 = sum(ϵr1) / 2 
ΛB1 = Ψ.(ϵ1 * eigvals(B1))
ΛD1 = Ψ.(ϵ1 * eigvals(D1))

# Now the eigenvalue condition satisfied. To show that the desired cluster synchronization is achieved let us solve the
# network. 

nodes1 = [Lorenz(σ=10, β=8/3, ρ=35) for i in 1 : n1] 
net1 = Network(nodes1, -ϵ1 * Ξ1, P)
sol1 = solvenet(net1, (0., 1000.))
t1, x1 = sol1.t, sol1.u 
d = dimension(net1.nodes[1])
plotres1 =  plot(layout=(4, 1)) 
plot!(t1, abs.(getwave(x1, 2, d) - getwave(x1, 5, 3)), label="e25", subplot=1)
plot!(t1, abs.(getwave(x1, 3, d) - getwave(x1, 4, 3)), label="e34", subplot=2)
plot!(t1, abs.(getwave(x1, 1, d) - getwave(x1, 2, 3)), label="e12", subplot=3)
plot!(t1, abs.(getwave(x1, 1, d) - getwave(x1, 3, 3)), label="e13", subplot=4)
display(plotres1)