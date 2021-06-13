using LightGraphs
using GraphPlot

n = 5   # Number of nodes in the network
graph = cycle_graph(n)
add_edge!(graph, 1, 3)
add_edge!(graph, 1, 4)
add_edge!(graph, 2, 5)
gplot(graph, nodelabel=1:nv(graph), layout = circular_layout)

W = collect(laplacian_matrix(graph)) .|> float

function scale(W0)
    W = copy(W0)
    for i in 1 : size(W, 1)
        W[i, :] ./=  W[i, i]
    end
    W
end
Ξ =  scale(W)

Q = [
    1 0 0 0 0;
    0 0 0 0 1;
    0 0 0 1 0;
    0 0 1 0 0;
    0 1 0 0 0
]

using LinearAlgebra
M = eigvecs(Q)

G = inv(M) * Q * M

B = Q[1 : 2, 1 : 2]

D = Q[3 : 5, 3 : 5]

λB =  eigvals(B)

λD = eigvals(D)

function sweepweight(w, W, M)
    λB = Vector{Vector{Float64}}(undef, length(w))
    λD = Vector{Vector{Float64}}(undef, length(w))
    for (i, wi) in enumerate(w)
        W[2, 5] = W[5, 2] =  -wi
        W[2, 2] -= sum(W[2, :])
        W[5, 5] -= sum(W[5, :])
        Ξ = scale(W)
        G = inv(M) * Ξ * M
        B = G[1 : 2, 1 : 2]
        D = G[3 : 5, 3 : 5]
        λB[i] = eigvals(B)
        λD[i] = eigvals(D)
    end
    λB, λD
end

w = collect(0 : 0.1 : 2)
λB, λD = sweepweight(w, W, M)

using Plots
theme(:default)
marker = (:circle, 3)
plot(w,  getindex.(λB, 1), marker=marker,  label="λ1tr")
plot!(w, getindex.(λB, 2), marker=marker,  label="λ2tr")
plot!(w, getindex.(λD, 1), marker=marker,  label="λ1syn")
plot!(w, getindex.(λD, 2), marker=marker,  label="λ2syn")
plot!(w, getindex.(λD, 3), marker=marker,  label="λ3syn")
vline!([1.], color="red", linestyle=:dash)
xlabel!("w")

using LyapunovExponents
ds = Lorenz(σ = 10, β = 8 / 3, ρ = 35)
nsteps = Int(3e4)   # Number of steps to calculate Lyapunov exponents
ntrsteps=Int(1e4)   # Number of steps transient steps before calcuting Lyapunov exponents
dt = 0.01           # Integration step size
lyaplorenz = lyapunovs(ds, nsteps = Int(3e4), ntrsteps=Int(1e4), dt=0.01)

P = [
    1 0 0;
    0 0 0;
    0 0 0
]

η = range(0, 100, length=21)
Ψ(ηi) = msf(ds, ηi, P)
Λ = Ψ.(η)
plot(η, Λ, marker=(:circle, 3))
hline!([0])
threshold = η[findfirst(Λ .< 0)]
@info "The critical value of η: $threshold"
vline!([threshold], label="Threshold: $threshold")

W = collect(laplacian_matrix(graph)) .|> float
wi = 1.5
W[2, 5] = W[5, 2] =  -wi
W[2, 2] -= sum(W[2, :])
W[5, 5] -= sum(W[5, :])
Ξ = scale(W)
G = inv(M) * Ξ * M
B = G[1 : 2, 1 : 2]
D = G[3 : 5, 3 : 5]
λB = eigvals(B)
λD = eigvals(D)
λ1tr = λB[1]
λ2syn = λD[2]
ϵ1 = threshold / λ1tr |> display
ϵ2 = threshold / λ2syn |> display

ϵ = 9.6
nodes = [Lorenz(σ = 10, β = 8 / 3, ρ=35) for i in 1 : n]
net = Network(nodes, -ϵ * Ξ, P)
sol = solvenet(net, (0., 1000.))
t, x = sol.t, sol.u

d = dimension(nodes[1])
getwave(x, i) = getindex.(x, (i - 1) * d + 1)

plt = plot(layout = (4, 1))
plot!(t, abs.(getwave(x, 2) -  getwave(x, 5)), label="e25", subplot=1)
plot!(t, abs.(getwave(x, 3) -  getwave(x, 4)), label="e34", subplot=2)
plot!(t, abs.(getwave(x, 1) -  getwave(x, 2)), label="e12", subplot=3)
plot!(t, abs.(getwave(x, 1) -  getwave(x, 3)), label="e13", subplot=4)
plt

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

