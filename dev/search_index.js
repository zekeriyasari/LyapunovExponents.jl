var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = LyapunovExponents","category":"page"},{"location":"#LyapunovExponents","page":"Home","title":"LyapunovExponents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for LyapunovExponents.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [LyapunovExponents]","category":"page"},{"location":"#LyapunovExponents.LyapunovExponents","page":"Home","title":"LyapunovExponents.LyapunovExponents","text":"A module to calculate lyapunov exponents.\n\n\n\n\n\n","category":"module"},{"location":"#LyapunovExponents.dimension-Tuple{LyapunovExponents.Dynamics}","page":"Home","title":"LyapunovExponents.dimension","text":"dimension(ds)\n\n\nReturns the dimension of ds \n\n\n\n\n\n","category":"method"},{"location":"#LyapunovExponents.lyapunovs-Tuple{LyapunovExponents.Dynamics}","page":"Home","title":"LyapunovExponents.lyapunovs","text":"lyapunovs(ds; kwargs...)\n\n\nComputes lyapunov exponents of ds.kwargs may include \n\nnsteps::Int = 3e4    Number of steps to calculate lypaunov exponents \nntrsteps::Int = 0   Number of transient steps for transient steps. Transient steps are taken before calculation of lypaunov exponents \ndt::Real = 0.01   Step size of integrator. \n\n\n\n\n\n","category":"method"},{"location":"#LyapunovExponents.msf-Tuple{LyapunovExponents.Dynamics, Real, AbstractMatrix{T} where T}","page":"Home","title":"LyapunovExponents.msf","text":"msf(ds, λ, P; kwargs...)\n\n\nMaster stability function (msf). Returns maximum lyapunov exponent. kwargs may include \n\nnsteps::Int = 3e4    Number of steps to calculate lypaunov exponents \nntrsteps::Int = 0   Number of transient steps for transient steps. Transient steps are taken before calculation of lypaunov exponents \ndt::Real = 0.01   Step size of integrator. \n\n\n\n\n\n","category":"method"}]
}
