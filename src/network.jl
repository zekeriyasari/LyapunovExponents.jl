# This file includes dynamical networks

export Network, solvenet, isstable, eigenmodes, unstablemodes, stablemodes

struct Network{T1<:Dynamics, T2<:AbstractMatrix, T3<:AbstractMatrix}
    nodes::Vector{T1}
    E::T2 
    P::T3 
end 
Network(::Type{T}, E::AbstractMatrix, P::AbstractMatrix) where {T<:Dynamics} = Network([T() for i in 1 : size(E, 1)], E, P)

function (net::Network)(dx, x, u, t) 
    for (k, idx) in enumerate(Iterators.partition(1 : length(x), size(net.P, 1)))
        net.nodes[k](view(dx, idx), view(x, idx), u, t) 
    end 
    dx .+= kron(net.E, net.P) * x 
end 

function solvenet(net::Network, tspan::Tuple, alg=Tsit5()) 
    prob = ODEProblem(net, vcat([node.x0 for node in net.nodes]...), tspan)
    solve(prob, SOLVER)
end

function eigenmodes(net::Network) 
    map(enumerate(eigvals(-net.E))) do (k, λi)
        λi => msf(net.nodes[k], λi, net.P)
    end
end 

function isstable(net::Network) 
    modepairs = eigenmodes(net) 
    _ = popfirst!(modepairs)
    all(getfield.(modepairs, :second) .< 0)
end 

unstablemodes(net::Network) = filter(modepair -> modepair.second > 0, eigenmodes(net))

stablemodes(net::Network) = filter(modepair -> modepair.second < 0, eigenmodes(net))
