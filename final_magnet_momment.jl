using LinearAlgebra
using Plots

# GEOMETRY
n_nodes = 11
nodes = [[Float64(i-1) * 0.1, 0.0] for i in 1:n_nodes]
edges = [[i, i+1] for i in 1:n_nodes-1]

# Material properties
k_spring = 1e4
damping = 0.95
B_field = [1.0, 1.0]
mag_moments = [10,20, 30, 50, 70 ,90 , 100]
fixed_nodes = [1]

# Simulation parameters
dt = 0.001
n_steps = 1000

# Function to run simulation for a given magnetic moment
function run_simulation(mag_moment, n_steps, dt)
    pos = copy(nodes)
    vel = [zeros(2) for _ in 1:n_nodes]
    force = [zeros(2) for _ in 1:n_nodes]
    mass = ones(n_nodes)
    pos_history = [copy(pos)]
    
    function compute_forces!(force, pos, mag_moment)
        fill!(force, [0.0, 0.0])
        
        for edge in edges
            i, j = edge
            r = pos[j] - pos[i]
            L = norm(r)
            L0 = norm(nodes[j] - nodes[i])
            if L > 1e-10
                F = k_spring * (L - L0) / L * r
                force[i] += F
                force[j] -= F
            end
        end
        
        force[end] += mag_moment * B_field
        
        for i in 1:n_nodes
            force[i] += [0.0, -9.81 * mass[i]]
        end
    end
    
    for step in 1:n_steps
        compute_forces!(force, pos, mag_moment)
        
        for i in 1:n_nodes
            if i ∉ fixed_nodes
                vel[i] += (force[i] / mass[i]) * dt
                vel[i] *= damping
                pos[i] += vel[i] * dt
            end
        end
        
        push!(pos_history, copy(pos))
    end
    
    return pos_history
end

# Run simulations for all magnetic moments
all_pos_histories = []
for mag_moment in mag_moments
    println("Running simulation for mag_moment = $mag_moment")
    pos_hist = run_simulation(mag_moment, n_steps, dt)
    push!(all_pos_histories, pos_hist)
end

# Create combined plot with all magnetic moment intensities
final_step = n_steps + 1

p = plot(xlim=(-0.1, 1.2), ylim=(-1.2, 0.5), 
         aspect_ratio=:equal, legend=:topright,
         title="Beam Deflection for Different Magnetic Moment Intensities",
         xlabel="x (m)", ylabel="y (m)",
         size=(1200, 800))

# Plot original beam
x_orig = [n[1] for n in nodes]
y_orig = [n[2] for n in nodes]
plot!(p, x_orig, y_orig, color=:black, linestyle=:dash, linewidth=2, 
      label="Original Position", alpha=0.5)

# Color palette for different magnetic moments
colors = palette(:tab10, length(mag_moments))

# Plot deformed beams for each magnetic moment
for (idx, mag_moment) in enumerate(mag_moments)
    pos_hist = all_pos_histories[idx]
    final_pos = pos_hist[final_step]
    
    x_def = [pt[1] for pt in final_pos]
    y_def = [pt[2] for pt in final_pos]
    
    plot!(p, x_def, y_def, color=colors[idx], linewidth=2.5, 
          label="μ = $mag_moment", alpha=0.8)
    scatter!(p, x_def, y_def, color=colors[idx], markersize=3, alpha=0.6, 
             label="")
end

# Highlight fixed node
scatter!(p, [nodes[1][1]], [nodes[1][2]], 
         color=:green, markersize=10, markershape=:square, 
         label="Fixed Node", legend=:bottomright)

savefig(p, "beam_deflection_all_moments.png")
println("Saved: beam_deflection_all_moments.png")
display(p)
