# In this file network simulation with changing network topology is performed. 

using LyapunovExponents 
using DifferentialEquations
using LightGraphs, GraphPlot 
using Plots; theme(:default) 

# Some helper functions 
plotgraph(graph::AbstractGraph) = gplot(graph, nodelabel=1:nv(graph), layout=circular_layout)
plotgraph(net::Network) = net.E |> SimpleGraph |> plotgraph 

# Construct a graph
n = 2       # Number of nodes 
ϵ = 20.     # Coupling strengths 
graph = cycle_graph(n)
Ξ = collect(laplacian_matrix(graph))
P = diagm([1, 0, 0])
net = Network(Lorenz, -ϵ * Ξ, P)
plotgraph(net) 

# Construct a callback 
actiontimes = [25., 50.] 
condition_delete(x, t, integ) = t ∈ actiontimes[1]  
condition_add(x, t, integ) = t ∈ actiontimes[2]  
action_delete!(integ) = deletelink!(integ.sol.prob.f.f, 1, 2) 
action_add!(integ) = addlink!(integ.sol.prob.f.f, 1, 2, ϵ) 
cb = CallbackSet(
    DiscreteCallback(condition_delete, action_delete!),
    DiscreteCallback(condition_add, action_add!)
) 

# Solve the network 
sol = solvenet(net, (0., 100.), callback=cb, tstops=actiontimes)
t, x = sol.t, sol.u 

# Plot the results 
d = dimension(net.nodes[1]) 
getwave(x, i) = getindex.(x, (i - 1) * d + 1)
plot(t, abs.(getwave(x, 1) - getwave(x, 2)), label="e12")
vline!(actiontimes, linestyle=:dash, linewidth=2, labels=["t", "s"])
