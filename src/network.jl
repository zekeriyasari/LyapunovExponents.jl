# This file includes dynamical networks

export Network, solvenet, isstable, eigenmodes, unstablemodes, stablemodes

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
