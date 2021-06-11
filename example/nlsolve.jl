# This file includes an example of nlsolve root finding. 

using NLsolve

# f(x) = (x - 1)^2 - 1
function f!(dx, x)
    dx[1] = (x[1] - 1)^2 - 1
end

sol = nlsolve(f!, [3.])
