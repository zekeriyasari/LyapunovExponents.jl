# # Updating Network Topology 
# In this notebook, we will consider mutating the network topology. For this purpose, let us construct a simple network
# and modify its topology step by step. 

# Load pacakges 
using LyapunovExponents 
using LightGraphs 
using GraphPlot 

# Let us define a helper fucntion to plot the network topology. 
plotgraph(graph::AbstractGraph) = gplot(graph, nodelabel=1:nv(graph), layout=circular_layout)
plotgraph(net::Network) = plotgraph(SimpleGraph(net.E))

# Construct a cyclic graph. 
n = 5   # Number of nodes 
graph = cycle_graph(n) 
plotgraph(graph)

# Construct a network from the graph. 
E = collect(laplacian_matrix(graph))
net = Network(Lorenz, E, diagm(ones(3))) 
plotgraph(net)

# Add and edge between the nodes 1 and 3 
updatelink!(net, 1, 3) 
plotgraph(net)

# Add and edge between the nodes 2 and 3 
updatelink!(net, 2, 5) 
plotgraph(net)

# Delete a connection between the nodes 2 and 3 
updatelink!(net, 2, 3, 0.)
plotgraph(net)