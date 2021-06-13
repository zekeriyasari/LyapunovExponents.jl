# This file includes dynamical networks

export Network, solvenet, isstable, eigenmodes, unstablemodes, stablemodes, updatelink!, scale_connection_matrix, 
    scale_connection_matrix!, haslink, addlink!, deletelink!, rewirelink!

"""
    $TYPEDEF 

Network of dynamical system nodes. 

# Fields 
    $TYPEDFIELDS 
"""
struct Network{T1<:Dynamics, T2<:AbstractMatrix, T3<:AbstractMatrix}
    "Dynamical system nodes"
    nodes::Vector{T1}
    "Outer coupling matrix"
    E::T2 
    "Inner coupling matrix"
    P::T3 
end 

Network(::Type{T}, E::AbstractMatrix, P::AbstractMatrix) where {T<:Dynamics} = Network([T() for i in 1 : size(E, 1)], E, P)

function (net::Network)(dx, x, u, t) 
    for (k, idx) in enumerate(Iterators.partition(1 : length(x), size(net.P, 1)))
        net.nodes[k](view(dx, idx), view(x, idx), u, t) 
    end 
    dx .+= kron(net.E, net.P) * x 
end 

"""
    $SIGNATURES 

Solves the network differential equation 
"""
function solvenet(net::Network, tspan::Tuple, alg=Tsit5(); solkwargs...) 
    prob = ODEProblem(net, vcat([node.x0 for node in net.nodes]...), tspan)
    solve(prob, SOLVER; solkwargs...)
end

"""
    $SIGNATURES

Returns eigenmodes, i.e. a vector of pairs of eigenvalues of the outer coupling matrix of `net` and corresponding maximum
Lyapunov exponent. 
"""
function eigenmodes(net::Network) 
    map(enumerate(eigvals(-net.E))) do (k, λi)
        λi => msf(net.nodes[k], λi, net.P)
    end
end 


"""
    $SIGNATURES

Returns true if `net` is stable. Stability analyis is performed via master stability function.
"""
function isstable(net::Network) 
    modepairs = eigenmodes(net) 
    _ = popfirst!(modepairs)
    all(getfield.(modepairs, :second) .< 0)
end 

"""
    $SIGNATURES

Returns unstable eigen modes of `net`. 
"""
unstablemodes(net::Network) = filter(modepair -> modepair.second > 0, eigenmodes(net))

"""
    $SIGNATURES

Returns stable eigen modes of `net`. 
"""
stablemodes(net::Network) = filter(modepair -> modepair.second < 0, eigenmodes(net))


# ------------------------------------- Topological modifications ----------------------------------------------------#


function islaplacian(E::AbstractMatrix)  
    ondiagidx = diagind(E)  # On-diagonal indices 
    offdiagidx = setdiff(1 : length(E), ondiagidx)  # Off-diagonal indices 
    all(E[ondiagidx] .> 0) && all(E[offdiagidx] .≤ 0)
end 

"""
    $SIGNATURES 

Returns true if `net` is constructed using laplacian matrix 
"""
islaplacian(net::Network) = islaplacian(net.E) 

"""
    $SIGNATURES 

Returns true if `net` is constructed using adjacency matrix 
"""
isadjacency(net::Network) = islaplacian(-net.E) 

"""
    $SIGNATURES 

Addes a link between the nodes `i` and `j` of the network `net`. `w` is the weight of the link to be updated.
"""
function updatelink!(net::Network, i::Int, j::Int, w::Real=islaplacian(net) ? -1 : 1) 
    islaplacian(net) && w > 0 && error("The network is laplacina matrix. Expected w to be nonnegative")
    isadjacency(net) && w < 0 && error("The network is adjacency matrix. Expected w to be nonpositive")
    issymmetric(net.E) || error("Expexted a symmetric matrix. Got, $(net.E)") 
    i == j && error("Expected different i and j, got i=$i, j = $j")
    E = net.E 
    colrange = 1 : size(E, 1) 
    E[i, j] = E[j, i] = w 
    E[i, i] = -sum(E[i, [k for k in colrange if k ≠ i]]) 
    E[j, j] = -sum(E[j, [k for k in colrange if k ≠ j]]) 
    net 
end 

"""
    $SIGNATURES

Adds a link between the nodes `i` and `j` with weight `w`. 
"""
addlink!(net::Network, i::Int, j::Int, w::Real=islaplacian(net) ? -1 : 1) = updatelink!(net, i, j, w) 

"""
    $SIGNATURES 

Deletes the links between the nodes `i` and `j` 
"""
deletelink!(net::Network, i::Int, j::Int) = updatelink!(net, i, j, 0.) 

""" 
    $SIGNATURES 

Removes the link from between the nodes `i` and `j` and wires the link between the nodes `k` and `l`. 
"""
function rewirelink!(net::Network, i::Int, j::Int, k::Int, l::Int)
    haslink(net, i, j) || error("The network does not have a link between $i and $j") 
    haslink(net, k, l) && error("The network has already a link between $l and $k") 
    updatelink!(net, k, l, net.E[i, j])  
end 

"""
    $SIGNATURES 

Returns true if there exists a link between the nodes `i` and `j`.
"""
haslink(net::E, i::Int, j::Int) 

"""
    $SIGNATURES 

Normalize the connection matrix of `net`. The connectin matrix is rescaled such that on-diagonal elements becones 1. 
"""
function scale_connection_matrix!(W0::AbstractMatrix) 
    for i in 1 : size(W0, 1)
        W0[i, :] ./=  W0[i, i]
    end
    W0
end 

scale_connection_matrix(W0::AbstractMatrix) = scale_connection_matrix!(copy(W0))

scale_connection_matrix!(net::Network) = scale_connection_matrix!(net.E) 
scale_connection_matrix(net::Network)  = scale_connection_matrix(net.E) 
