using ForwardDiff 
using LyapunovExponents 

ds = Lorenz() 
d = dimension(ds)
J0 = zeros(d, d)
x0 = ds.x0
J = ds(J0, x0, nothing, 0.)

function defjacob(ds::Dynamics)
    (ds::Dynamics)(J::AbstractMatrix, x::AbstractVector, u, t) = 
        ForwardDiff.jacobian!(J, (dx, x) -> ds(dx, x, nothing, 0.), deepcopy(ds.x0), ds.x0)
end 

