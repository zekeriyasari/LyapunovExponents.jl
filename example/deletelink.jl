# # Topological Control of Dynamical Networks 

# In this notebook, we will consider topological control of dynamical networks. 

# First load required packages and  define some helper functions 

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
plotgraph(graph) = gplot(graph, nodelabel=1:nv(graph), layout=circular_layout)
getwave(x, i, d) = getindex.(x, (i - 1) * d + 1) 

# ## Master Stability Function 

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

# ## Network 1

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

# ## Network 2 - Delete a Link 

# Note that the eigeen value condtion is not satisfied since the first eigenvalue of the transversal subspace and least
# nontrivial eigenvalue of the synchronization subspace overlaps. To make them seperated, let us delete the link between
# the nodes 1 and 3. When the link changes, the symmetry of the network also changes.

graph2 = copy(graph1) 
rem_edge!(graph2, 1, 3) 
plotgraph(graph2) |> display 
Ξ2 = collect(laplacian_matrix(graph2))
Q2 = [
    0 0 0 0 1; 
    0 0 0 1 0; 
    0 0 1 0 0; 
    0 1 0 0 0;
    1 0 0 0 0
]
M2 = eigvecs(Q2)
G2 = inv(M2) * Ξ2 * M2 
B2 = G2[1 : 2, 1 : 2] 
D2 = G2[3 : 5, 3 : 5] 
ploteigs2 = plot() 
λB2 = eigvals(B2) 
λD2 = eigvals(D2) 
scatter!(λB2, zeros(length(λB2)), marker=mt, label="tr", xlabel="λ")
scatter!(λD2, zeros(length(λD2)), marker=ms, label="syn", xlabel="λ")
display(ploteigs2)

ϵr2  = ηtr ./  [eigvals(B2)[1], eigvals(D2)[2]] 
ϵ2 = sum(ϵr2) / 2 
ΛB2 = Ψ.(ϵ2 * eigvals(B2))
ΛD2 = Ψ.(ϵ2 * eigvals(D2))

# Now the eigenvalue condition satisfied. To show that the desired cluster synchronization is achieved let us solve the
# network. 

nodes2 = [Lorenz(σ=10, β=8/3, ρ=35) for i in 1 : n1] 
net2 = Network(nodes2, -ϵ2 * Ξ2, P)
sol2 = solvenet(net2, (0., 1000.))
t2, x2 = sol2.t, sol2.u 
d = dimension(net2.nodes[1])
plotres2 =  plot(layout=(4, 1)) 
plot!(t2, abs.(getwave(x2, 1, d) - getwave(x2, 5, d)), label="e15", subplot=1)
plot!(t2, abs.(getwave(x2, 2, d) - getwave(x2, 4, d)), label="e24", subplot=2)
plot!(t2, abs.(getwave(x2, 1, d) - getwave(x2, 3, d)), label="e13", subplot=3)
plot!(t2, abs.(getwave(x2, 2, d) - getwave(x2, 3, d)), label="e23", subplot=4)
display(plotres2)

# Check synchronization stability of the previous network 

ΛB2 = Ψ.(ϵ2 * eigvals(B1))
ΛD2 = Ψ.(ϵ2 * eigvals(D1))

# Note that the network before the removal of the link between the node 1 and 3 is fully synchronized. Let see this
# case, 

nodes1 = [Lorenz(σ=10, β=8/3, ρ=35) for i in 1 : n1] 
net1 = Network(nodes1, -ϵ2 * Ξ1, P)  # First network before the link is removed. 
sol1 = solvenet(net1, (0., 1000.))
t1, x1 = sol1.t, sol1.u 
d = dimension(net1.nodes[1])
plotres1 =  plot(layout=(4, 1)) 
plot!(t1, abs.(getwave(x1, 1, d) - getwave(x1, 2, d)), label="e12", subplot=1)
plot!(t1, abs.(getwave(x1, 1, d) - getwave(x1, 2, d)), label="e13", subplot=2)
plot!(t1, abs.(getwave(x1, 1, d) - getwave(x1, 3, d)), label="e14", subplot=3)
plot!(t1, abs.(getwave(x1, 1, d) - getwave(x1, 4, d)), label="e15", subplot=4)
display(plotres1)

# Now let us see how cluster pattern of the network changes when the link is removed 
nodes3 = [Lorenz(σ=10, β=8/3, ρ=35) for i in 1 : n1] 
net3 = Network(nodes3, -ϵ2 * Ξ1, P)  # First network before the link is removed. 

actiontimes = [100, 500.] 
condition_delete(x, t, integ) = t ∈ actiontimes[1]  
condition_add(x, t, integ) = t ∈ actiontimes[2]  
action_delete!(integ) = deletelink!(integ.sol.prob.f.f, 1, 3) 
action_add!(integ) = addlink!(integ.sol.prob.f.f, 1, 3, ϵ2) 
cb = CallbackSet(
    DiscreteCallback(condition_delete, action_delete!),
    DiscreteCallback(condition_add, action_add!)
) 

sol3 = solvenet(net3, (0., 1000.), callback=cb, tstops=actiontimes)
t3, x3 = sol3.t, sol3.u 
d = dimension(net3.nodes[1])
plotres3 =  plot(layout=(4, 1)) 
plot!(t3, abs.(getwave(x3, 1, d) - getwave(x3, 5, d)), label="e15", subplot=1)
plot!(t3, abs.(getwave(x3, 2, d) - getwave(x3, 4, d)), label="e24", subplot=2)
plot!(t3, abs.(getwave(x3, 1, d) - getwave(x3, 3, d)), label="e13", subplot=3)
plot!(t3, abs.(getwave(x3, 2, d) - getwave(x3, 3, d)), label="e23", subplot=4)
for i in 1 : 4 
    vline!(actiontimes, linestyle=:dash, linewidth=2, subplot=i)
end 
display(plotres3)

# Note that the synchronization pattern of the network changes from full synchronization to cluster synchronization and
# full synchronization again. 

# ## Network 3 - Revire a Link 

# Now we can investigate the effect of deleting a link between a pair of nodes and rewiring the link between different
# pair of nodes. Consider the following network 

graph3 = copy(graph2) 
rem_edge!(graph3, 1, 4) 
add_edge!(graph3, 1, 3)
plotgraph(graph3) |> display 
Ξ3 = collect(laplacian_matrix(graph3))
Q3 = [
    0 1 0 0 0; 
    1 0 0 0 0; 
    0 0 0 0 1; 
    0 0 0 1 0;
    0 0 1 0 0
]
M3 = eigvecs(Q3)
G3 = inv(M3) * Ξ3 * M3 
B3 = G3[1 : 2, 1 : 2] 
D3 = G3[3 : 5, 3 : 5] 
ploteigs3 = plot() 
λB3 = eigvals(B3) 
λD3 = eigvals(D3) 
scatter!(λB3, zeros(length(λB3)), marker=mt, label="tr", xlabel="λ")
scatter!(λD3, zeros(length(λD3)), marker=ms, label="syn", xlabel="λ")
display(ploteigs3)

ϵr3  = ηtr ./  [eigvals(B3)[1], eigvals(D3)[2]] 
ϵ3 = sum(ϵr3) / 2 
ΛB3 = Ψ.(ϵ3 * eigvals(B3))
ΛD3 = Ψ.(ϵ3 * eigvals(D3))

# Check synchronization stability of the previous network 

# Now let us see how cluster pattern of the network changes when the link is removed 
nodes4 = [Lorenz(σ=10, β=8/3, ρ=35) for i in 1 : n1] 
net4 = Network(nodes4, -ϵ3 * Ξ2, P)  # First network before the link is removed. 

actiontimes = [200, 500.] 
condition_delete(x, t, integ) = t ∈ actiontimes[1]  
condition_add(x, t, integ) = t ∈ actiontimes[2]  
action_delete!(integ) = deletelink!(integ.sol.prob.f.f, 1, 4) 
action_add!(integ) = addlink!(integ.sol.prob.f.f, 1, 3, ϵ3) 
cb = CallbackSet(
    DiscreteCallback(condition_delete, action_delete!),
    DiscreteCallback(condition_add, action_add!)
) 

sol4 = solvenet(net4, (0., 1000.), callback=cb, tstops=actiontimes)
t4, x4 = sol4.t, sol4.u 
d = dimension(net4.nodes[1])
plotres3 =  plot(layout=(6, 1)) 
plot!(t4, abs.(getwave(x4, 1, d) - getwave(x4, 5, d)), label="e15", subplot=1)
plot!(t4, abs.(getwave(x4, 2, d) - getwave(x4, 4, d)), label="e24", subplot=2)
plot!(t4, abs.(getwave(x4, 1, d) - getwave(x4, 2, d)), label="e12", subplot=3)
plot!(t4, abs.(getwave(x4, 3, d) - getwave(x4, 5, d)), label="e35", subplot=4)
plot!(t4, abs.(getwave(x4, 1, d) - getwave(x4, 3, d)), label="e13", subplot=5)
plot!(t4, abs.(getwave(x4, 2, d) - getwave(x4, 3, d)), label="e23", subplot=6)
for i in 1 : 4 
    vline!(actiontimes, linestyle=:dash, linewidth=2, subplot=i)
end 
display(plotres3)
