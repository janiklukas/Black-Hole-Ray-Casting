"""
The function `reduce_angles` takes in a state vector `X::MVector{8,Float64}` and returns a new state with the angular coordinates reset to the correct domains
"""
function reduce_angles(X::MVector{8,Float64})
    if X[3] < 0
        X[3] = abs(X[3])
        X[4] += π
    end

    X[3] = mod(X[3], 2π)

    if X[3] > π
        X[3] = 2π - X[3]
        X[4] += π
    end
    
    X[4] = mod(X[4], 2π)

    return X
end

"""
The function `Schwarzschild` takes in a state vector `X::MVector{8,Float64}` and returns the right-hand side of the geodesic ODE system
"""
function Schwarzschild(X::MVector{8,Float64})
    X = reduce_angles(X)
    dX = MVector{8,Float64}(zeros(8))

    dX[1] = X[5]; dX[2] = X[6]; dX[3] = X[7]; dX[4] = X[8]
    dX[5] = 1/(X[2]-X[2]^2) * X[5] * X[6]
    #dX[6] = (1-X[2])/(2*X[2]^3) * X[5]^2 + 1/(2*X[2]^2-2*X[2]) * X[6]^2 + (X[2]-1) * (X[7]^2 + sin(X[3])^2 * X[8]^2) #without lightlike condition
    dX[6] = (2*X[2]-3)/2 * (X[7]^2 + sin(X[3])^2 * X[8]^2) #with lightlike condition
    dX[7] = -2/X[2] * X[6] * X[7] + sin(X[3])*cos(X[3]) * X[8]^2
    dX[8] = -2/X[2] * X[6] * X[8] - 2*cot(X[3]) * X[7] * X[8] #watch out for cotangent NaNs!

    return dX
end

"""
The function `step_euler` takes in a state vector `X::MVector{8,Float64}`, advances it by the time step `Δt::Float64` using the Euler method and returns the new state
"""
function step_euler(X::MVector{8,Float64}, Δt::Float64)
    k1 = Schwarzschild(X)
    return X + Δt * k1
end

"""
The function `step_RK2` takes in a state vector `X::MVector{8,Float64}`, advances it by the time step `Δt::Float64` using the 2nd order Runge-Kutta method and returns the new state
"""
function step_RK2(X::MVector{8,Float64}, Δt::Float64)
    k1 = Schwarzschild(X)
    k2 = Schwarzschild(X + Δt * k1)
    return X + Δt/2 * (k1 + k2)
end

"""
The function `step_RK4` takes in a state vector `X::MVector{8,Float64}`, advances it by the time step `Δt::Float64` using the 4th order Runge-Kutta method and returns the new state
"""
function step_RK4(X::MVector{8,Float64}, Δt::Float64)
    k1 = Schwarzschild(X)
    k2 = Schwarzschild(X + Δt/2 * k1)
    k3 = Schwarzschild(X + Δt/2 * k2)
    k4 = Schwarzschild(X + Δt * k3)
    return X + Δt/6 * (k1 + 2*(k2 + k3) + k4)
end

"""
The function `step_state` selects an ODE solver based on the global variable `solver::String`
"""
function step_state(X::MVector{8,Float64}, Δt::Float64)
    if solver == "euler"
        return step_euler(X, Δt)
    elseif solver == "RK2"
        return step_RK2(X, Δt)
    elseif solver == "RK4"
        return step_RK4(X, Δt)
    else
        @error "`solver` must be either `euler`, `RK2` or `RK4`"
    end
end

"""
The function `update_ray!` advances a ray matrix `ray::MMatrix{2,8,Float64}` by the time step `Δt::Float64`
"""
function update_ray!(ray::MMatrix{2,8,Float64}, Δt::Float64)
    ray[1, :] = reduce_angles(ray[2, :])
    ray[2, :] = reduce_angles(step_state(ray[2, :], Δt))
end

"""
The function `sign_flip` tests if the θ = π/2 plane is crossed between two angles `θ_prev::Float64` and `θ_curr::Float64`
"""
function sign_flip(θ_prev::Float64, θ_curr::Float64)
    if (θ_prev - π/2) * (θ_curr - π/2) < 0
        return true
    else
        return false
    end
end

"""
The function `detect_hit` checks if a given ray that crosses the θ = π/2 plane hits the accretion disc using an interval search
### Arguments:
- `ray::MMatrix{2,8,Float64}`: the light ray matrix
- `ρ::Float64`: the accretion disc radius
- `Δt::Float64`: the time step size
- `iters::Int64`: the number of interval search iterations
"""
function detect_hit(ray::MMatrix{2,8,Float64}, Δt::Float64, ρ::Float64, iters::Int64)
    #Refine the trajectory using an interval search
    X_prev = ray[1, :]
    X_curr = ray[2, :]
    for i in 1:iters
        τ = Δt / (2^i)
        X_mid = step_state(X_prev, τ)
        
        if sign_flip(X_prev[3], X_mid[3])
            X_curr = X_mid
        elseif sign_flip(X_mid[3], X_curr[3])
            X_prev = X_mid
        else
            @warn "hit detection failed"
        end
    end
    
    #Check for disc collision
    if X_curr[2] <= 1 + ρ
        return true
    else
        return false
    end
end

"""
The function `ray_cast` evolves a given initial ray, checking for hits after each step and returning 1 if the disc is hit and 0 otherwise, with a chance to return 2 for a background star
### Arguments:
- `ray::MMatrix{2,8,Float64}`: the initial light ray
- `Δt::Float64`: the time step size
- `t_max::Float64`: the time evolution cutoff
- `ρ::Float64`: the accretion disc radius
- `r_max::Float64`: the cutoff for the radial coordinate
- `iters::Int64`: the number of iterations for the hit detection function
- `star_prob::Float64`: the probability for adding background stars (optional)
"""
function ray_cast(ray::MMatrix{2,8,Float64}, Δt::Float64, t_max::Float64, ρ::Float64, r_max::Float64, iters::Int64; star_prob::Float64=0.0)
    steps = ceil(t_max / Δt)

    for _ in 1:steps
        #Check for horizon collisions
        if ray[2, 2] < 1 || abs(ray[2, 2] - 1) < 1e-3
            return 0
        end
        
        #Check for boundary exit and excessive velocity
        if abs(ray[2, 2]) > r_max || abs(ray[2, 6]) > 10
            break
        end

        #Check for disc collisions
        if sign_flip(ray[1, 3], ray[2, 3])
            if detect_hit(ray, Δt, ρ, iters)
                return 1
            end
        end

        update_ray!(ray, Δt)
    end

    #Randomly add some background stars
    if rand() < star_prob
        return 2
    else
        return 0
    end
end