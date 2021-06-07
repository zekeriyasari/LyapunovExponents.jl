# This file inludes the master stability function 

msf(ds::Dynamics, E::AbstractMatrix, P::AbstractMatrix) = map(λ -> _msf(ds, λ, P), eigvals(E))

function _msf(ds::Dynamics, λ::Real, P::AbstractMatrix) 
    d = dimension(ds) 
    J0 = ds(zeros(d, d), ds.x0, nothing, ds.t) 
end

function augment(ds::Dynamics, J::AbstractMatrix, λ::Real, P::AbstractMatrix)
    df = augment(ds, J) 
    function f(dΦ, Φ, J, t) 
        df(dΦ, Φ, J, t)
        Δ = @view Φ[:, 2 : end] 
        dΔ = @view dΦ[:, 2 : end]
        dΔ .-=  λ * P * Δ  # This part comes from coupling part.
    end 
end 
