"""
The function `generate_image` initializes the light rays, evolves them and shows the final result
### Arguments:
- `camera::Camera`: the camera parameter object
- `px::Int64`: the number of pixels in each dimension
- `Δt::Float64`: the time step size
- `t_max::Float64`: the time evolution cutoff
- `ρ::Float64`: the accretion disc radius
- `r_max::Float64`: the cutoff for the radial coordinate
- `iters::Int64`: the number of iterations for the hit detection function
- `size::Int64`: the side length of the output image
- `star_prob::Float64`: the probability for adding background stars (optional)
"""
function generate_image(camera::Camera, px::Int64, Δt::Float64, t_max::Float64, ρ::Float64, r_max::Float64, iters::Int64, size::Int64; star_prob::Float64=0.0)
    rays = rays_Schwarzschild(camera, px)
    img = map(copy, fill(RGB(0, 0, 0), (px,px)))

    for i in 1:px
        for j in 1:px
            source = ray_cast(rays[i, j], Δt, t_max, ρ, r_max, iters; star_prob)
            if source == 1
                img[i, j] = RGB(1, 0.3, 0)
            elseif source == 2
                img[i, j] = RGB(1, 1, 1)
            end
        end
    end

    img = imresize(img, (size,size))
    return RGB.(img)
end

"""
The function `animate_timesteps` produces an animation with varying time resolution
### Arguments:
- `camera::Camera`: the camera parameter object
- `px::Int64`: the number of pixels in each dimension
- `Δt_array::Float64`: the array of time step sizes
- `t_max::Float64`: the time evolution cutoff
- `ρ::Float64`: the accretion disc radius
- `r_max::Float64`: the cutoff for the radial coordinate
- `iters::Int64`: the number of iterations for the hit detection function
- `size::Int64`: the side length of the output image
- `fps::Int64`: the animation frames per second
- `gif_length::Float64`: the animation length
"""
function animate_timesteps(camera::Camera, px::Int64, Δt_array::Vector{Float64}, t_max::Float64, ρ::Float64, r_max::Float64, iters::Int64, size::Int64, fps::Int64, gif_length::Float64; star_prob::Float64=0.0)
    rays = rays_Schwarzschild(camera, px)

    anim = @animate for Δt in Δt_array
        plot(generate_image(camera, px, Δt, t_max, ρ, r_max, iters, size; star_prob))
        plot!(title = "Δt = " * string(round(Δt, sigdigits = 3)))
    end every round(Int, 1/(gif_length*fps) * length(Δt_array))
    
    display(gif(anim, "Animation.gif", fps = fps))
end