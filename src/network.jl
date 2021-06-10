# This file includes dynamical networks

export Network, solvenet, isstable, eigenmodes

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

# FIXME:The number of eigenmodes does not mathc the number of eigenvalues.  
function eigenmodes(net::Network) 
    λ = sort!(eigvals(-net.E))
    @show λ
    d = Dict{eltype(λ), Float64}()
    for (k, λi) in enumerate(λ)
        d[λi] = msf(net.nodes[k], λi, net.P)
    end 
    d
end 

isstable(net::Network) = all(value(eigenmodes) .≤ 0)
