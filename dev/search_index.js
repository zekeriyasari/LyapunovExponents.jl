var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = LyapunovExponents","category":"page"},{"location":"#LyapunovExponents","page":"Home","title":"LyapunovExponents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for LyapunovExponents.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [LyapunovExponents]","category":"page"},{"location":"#LyapunovExponents.LyapunovExponents","page":"Home","title":"LyapunovExponents.LyapunovExponents","text":"A module to calculate lyapunov exponents.\n\n\n\n\n\n","category":"module"},{"location":"#LyapunovExponents.Network","page":"Home","title":"LyapunovExponents.Network","text":"struct Network{T1<:LyapunovExponents.Dynamics, T2<:(AbstractMatrix{T} where T), T3<:(AbstractMatrix{T} where T)}\n\nNetwork of dynamical system nodes. \n\nFields\n\nnodes::Vector{T1} where T1<:LyapunovExponents.Dynamics\nDynamical system nodes\nE::AbstractMatrix{T} where T\nOuter coupling matrix\nP::AbstractMatrix{T} where T\nInner coupling matrix\n\n\n\n\n\n","category":"type"},{"location":"#LyapunovExponents.addlink!","page":"Home","title":"LyapunovExponents.addlink!","text":"addlink!(net, i, j)\naddlink!(net, i, j, w)\n\n\nAdds a link between the nodes i and j with weight w. \n\n\n\n\n\n","category":"function"},{"location":"#LyapunovExponents.deletelink!-Tuple{Network, Int64, Int64}","page":"Home","title":"LyapunovExponents.deletelink!","text":"deletelink!(net, i, j)\n\n\nDeletes the links between the nodes i and j \n\n\n\n\n\n","category":"method"},{"location":"#LyapunovExponents.dimension-Tuple{LyapunovExponents.Dynamics}","page":"Home","title":"LyapunovExponents.dimension","text":"dimension(ds)\n\n\nReturns the dimension of ds \n\n\n\n\n\n","category":"method"},{"location":"#LyapunovExponents.eigenmodes-Tuple{Network}","page":"Home","title":"LyapunovExponents.eigenmodes","text":"eigenmodes(net)\n\n\nReturns eigenmodes, i.e. a vector of pairs of eigenvalues of the outer coupling matrix of net and corresponding maximum Lyapunov exponent. \n\n\n\n\n\n","category":"method"},{"location":"#LyapunovExponents.haslink-Tuple{Network, Int64, Int64}","page":"Home","title":"LyapunovExponents.haslink","text":"haslink(net, i, j)\n\n\nReturns true if there exists a link between the nodes i and j.\n\n\n\n\n\n","category":"method"},{"location":"#LyapunovExponents.isadjacency-Tuple{Network}","page":"Home","title":"LyapunovExponents.isadjacency","text":"isadjacency(net)\n\n\nReturns true if net is constructed using adjacency matrix \n\n\n\n\n\n","category":"method"},{"location":"#LyapunovExponents.islaplacian-Tuple{Network}","page":"Home","title":"LyapunovExponents.islaplacian","text":"islaplacian(net)\n\n\nReturns true if net is constructed using laplacian matrix \n\n\n\n\n\n","category":"method"},{"location":"#LyapunovExponents.isstable-Tuple{Network}","page":"Home","title":"LyapunovExponents.isstable","text":"isstable(net)\n\n\nReturns true if net is stable. Stability analyis is performed via master stability function.\n\n\n\n\n\n","category":"method"},{"location":"#LyapunovExponents.lyapunovs-Tuple{LyapunovExponents.Dynamics}","page":"Home","title":"LyapunovExponents.lyapunovs","text":"lyapunovs(ds; kwargs...)\n\n\nComputes lyapunov exponents of ds.kwargs may include \n\nnsteps::Int = 3e4    Number of steps to calculate lypaunov exponents \nntrsteps::Int = 0   Number of transient steps for transient steps. Transient steps are taken before calculation of lypaunov exponents \ndt::Real = 0.01   Step size of integrator. \n\n\n\n\n\n","category":"method"},{"location":"#LyapunovExponents.msf-Tuple{LyapunovExponents.Dynamics, Real, AbstractMatrix{T} where T}","page":"Home","title":"LyapunovExponents.msf","text":"msf(ds, η, P; kwargs...)\n\n\nMaster stability function (msf). Returns maximum lyapunov exponent of the network defined as \n\ndotx_i = f(x)_i + sum_j = 1^n epsilon_ij P x_j quad i = 1 2 ldots n\n\nwhere x_i in mathbbR^d is the state vector of the i^th node, f is the vector-valued function defining the individual node dynamics. epsilon_ij geq 0 if there is a connection between the nodes i and j. epsilon_ij = 0, if there is no connection between the nodes i and j. E = epsilon_ij determines he network topology. P = diag(p_1 p_2 ldots p_d) determines the variables by which the nodes are conected to each other. The network can be expressed more compaclty as\n\ndotX = F(X) + (E otimes P) X \n\nwhere X = x_1 x_2 ldots x_n F(X) = f(x_1) f_(x_2) ldots f_(x_n) E = epsilon_ij P = diag(p_1 p_2 ldots p_n)\n\nkwargs are \n\nnsteps::Int = 3e4    Number of steps to calculate lypaunov exponents \nntrsteps::Int = 0   Number of transient steps for transient steps. Transient steps are taken before calculation of lypaunov exponents \ndt::Real = 0.01   Step size of integrator. \n\nnote: Note\nNote that epsilon_ij geq 0 j neq iEhas nonpositive eigenvalues Let\\etadenote the eigenvalue ofE``, than we have,     nu_n leq nu_n - 1 leq ldots leq nu_2  nu_1 = 0\n\n\n\n\n\n","category":"method"},{"location":"#LyapunovExponents.rewirelink!-Tuple{Network, Int64, Int64, Int64, Int64}","page":"Home","title":"LyapunovExponents.rewirelink!","text":"rewirelink!(net, i, j, k, l)\n\n\nRemoves the link from between the nodes i and j and wires the link between the nodes k and l. \n\n\n\n\n\n","category":"method"},{"location":"#LyapunovExponents.scale_connection_matrix!-Tuple{AbstractMatrix{T} where T}","page":"Home","title":"LyapunovExponents.scale_connection_matrix!","text":"scale_connection_matrix!(W0)\n\n\nNormalize the connection matrix of net. The connectin matrix is rescaled such that on-diagonal elements becones 1. \n\n\n\n\n\n","category":"method"},{"location":"#LyapunovExponents.solvenet","page":"Home","title":"LyapunovExponents.solvenet","text":"solvenet(net, tspan)\nsolvenet(net, tspan, alg; solkwargs...)\n\n\nSolves the network differential equation \n\n\n\n\n\n","category":"function"},{"location":"#LyapunovExponents.stablemodes-Tuple{Network}","page":"Home","title":"LyapunovExponents.stablemodes","text":"stablemodes(net)\n\n\nReturns stable eigen modes of net. \n\n\n\n\n\n","category":"method"},{"location":"#LyapunovExponents.unstablemodes-Tuple{Network}","page":"Home","title":"LyapunovExponents.unstablemodes","text":"unstablemodes(net)\n\n\nReturns unstable eigen modes of net. \n\n\n\n\n\n","category":"method"},{"location":"#LyapunovExponents.updatelink!","page":"Home","title":"LyapunovExponents.updatelink!","text":"updatelink!(net, i, j)\nupdatelink!(net, i, j, w)\n\n\nAddes a link between the nodes i and j of the network net. w is the weight of the link to be updated.\n\n\n\n\n\n","category":"function"}]
}
