"""
The `Camera` struct contains the relevant camera parameters, namely its coordinates `r::Float64` and `θ::Float64` as well as the distance to the screen `dist::Float64` and the field of view `fov::Float64`
"""
struct Camera
    r::Float64
    θ::Float64
    dist::Float64
    fov::Float64
end

"""
The function rays_Cartesian takes in a set `camera::Camera` of camera parameters and a pixel number `px::Int64` and returns an array of `px^2` initial rays in Cartesian screen coordinates
"""
function rays_Cartesian(camera::Camera, px::Int64)
    #Calculate screen size, pixel size and pixel number:
    s = 2 * camera.dist * tan(camera.fov/2)
    p = s / px

    #Find pixel locations
    x_coord = camera.r - camera.dist
    yz_coords = range(start = -s/2 + p/2, step = p, length = px)
    locs = [[x_coord, y, z] for y in yz_coords for z in yz_coords]

    #Find & normalize initial velocities
    vels = locs - fill([camera.r, 0, 0], (px^2,))
    for i in 1:px^2
        vels[i] = vels[i] / norm(vels[i])
    end
    
    #Assemble output (watch out for identical instances created by fill!)
    locs = reshape(locs, (px, px))
    vels = reshape(vels, (px, px))

    output = [MMatrix{2,8,Float64}(zeros((2, 8))) for _ in 1:px^2]
    output = reshape(output, (px, px))
    for i in 1:px
        for j in 1:px
            output[i, j][1, :] = [0; camera.r; camera.θ; 0; 0; vels[i, j]]
            output[i, j][2, :] = [0; locs[i, j]; 0; vels[i, j]]
        end
    end

    return output
end

"""
The function `convert_to_spherical` takes in a state vector `X::MVector{8,Float64}` in Cartesian screen coordinates and a `camera::Camera` object and transforms the spatial part of `X` to flat spherical coordinates
"""
function convert_to_spherical(X::MVector{8,Float64}, camera::Camera)
    loc = copy(X[2:4])
    vel = copy(X[6:8])

    #Rotate to standard Cartesian coordinates (about the y-axis):
    inc = camera.θ - π/2
    Rot = [cos(inc) 0 -sin(inc); 0 1 0; sin(inc) 0 cos(inc)]
    loc = Rot * loc
    vel = Rot * vel
    
    #Transform coordinates:
    r = norm(loc)
    θ = acos(loc[3]/r)
    ϕ = sign(loc[2]) * acos(loc[1] / norm(loc[1:2]))
    X[2:4] = [r, θ, ϕ]
    
    #Transform velocity:
    Jac = [sin(θ)*cos(ϕ) sin(θ)*sin(ϕ) cos(θ); r*cos(θ)*cos(ϕ) r*cos(θ)*sin(ϕ) -r*sin(θ); -r*sin(θ)*sin(ϕ) r*sin(θ)*cos(ϕ) 0]
    X[6:8] = Jac * vel

    #Fix angular domains:
    X = reduce_angles(X)

    return X
end

"""
#The function speed_Schwarzschild takes in a state vector `X::MVector{8,Float64}` and returns the spatial dot product of the velocity divided by g_00
"""
function speed_Schwarzschild(X::MVector{8,Float64})
    g_spatial = Diagonal([1/(1-1/X[2]), X[2]^2, X[2]^2*sin(X[3])^2])
    
    return dot(X[6:8], g_spatial, X[6:8])  / (1 - 1/X[2])
end

"""
The function rays_Schwarzschild produces an array of lightlike initial rays in Schwarzschild coordinates for a given `camera::Camera` object and pixel number `px::Int`
"""
function rays_Schwarzschild(camera::Camera, px::Int)
    rays = rays_Cartesian(camera, px)

    for i in 1:px
        for j in 1:px
            #Convert to spherical coordinates
            ray = rays[i, j]
            ray[1, :] = convert_to_spherical(ray[1, :], camera)
            ray[2, :] = convert_to_spherical(ray[2, :], camera)
            
            #Enforce lightlike condition
            ray[1, 5] = sqrt(speed_Schwarzschild(ray[1, :]))
            ray[2, 5] = sqrt(speed_Schwarzschild(ray[2, :]))
        end
    end

    return rays
end