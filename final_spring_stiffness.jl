using LinearAlgebra
using Plots

# GEOMETRY
n_nodes = 11
nodes = [[Float64(i-1) * 0.1, 0.0] for i in 1:n_nodes]
edges = [[i, i+1] for i in 1:n_nodes-1]

# Material properties
mag_moment = 50.0  # Fixed magnetic moment
damping = 0.9  # Fixed damping
B_field = [1.0, 1.0]
spring_stiffnesses = [1e3, 1e4, 1e5]  # 10^3, 10^4, 10^5
fixed_nodes = [1]

# Simulation parameters
dt = 0.001
n_steps = 1000

# Function to run simulation for a given spring stiffness
function run_simulation(k_spring, n_steps, dt)
    pos = copy(nodes)
    vel = [zeros(2) for _ in 1:n_nodes]
    force = [zeros(2) for _ in 1:n_nodes]
    mass = ones(n_nodes)
    pos_history = [copy(pos)]
    
    function compute_forces!(force, pos)
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
        compute_forces!(force, pos)
        
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

# Run simulations for all spring stiffnesses
all_pos_histories = []
for k_spring in spring_stiffnesses
    println("Running simulation for k_spring = $k_spring")
    pos_hist = run_simulation(k_spring, n_steps, dt)
    push!(all_pos_histories, pos_hist)
end

# Create single plot with all spring stiffnesses overlaid
final_step = n_steps + 1

p = plot(xlim=(-0.1, 1.2), ylim=(-1.2, 0.5), 
         aspect_ratio=:equal, legend=:topright,
         title="Beam Deflection for Different Spring Stiffnesses (μ = $(mag_moment), c = $damping)",
         xlabel="x (m)", ylabel="y (m)",
         size=(1000, 700))

# Plot original beam
x_orig = [n[1] for n in nodes]
y_orig = [n[2] for n in nodes]
plot!(p, x_orig, y_orig, color=:black, linestyle=:dash, linewidth=2.5, 
      label="Original Position", alpha=0.6)

# Color palette for different spring stiffnesses
colors = [:blue, :red, :green]

# Plot deformed beams for each spring stiffness
for (idx, k_spring) in enumerate(spring_stiffnesses)
    pos_hist = all_pos_histories[idx]
    final_pos = pos_hist[final_step]
    
    # Determine exponent for label
    exponent = Int(log10(k_spring))
    
    x_def = [pt[1] for pt in final_pos]
    y_def = [pt[2] for pt in final_pos]
    
    plot!(p, x_def, y_def, color=colors[idx], linewidth=3, 
          label="k = 10^$exponent N/m", alpha=0.8)
    scatter!(p, x_def, y_def, color=colors[idx], markersize=4, alpha=0.6, 
             label="")
end

# Highlight fixed node
scatter!(p, [nodes[1][1]], [nodes[1][2]], 
         color=:green, markersize=10, markershape=:square, 
         label="Fixed Node")

savefig(p, "beam_deflection_stiffness_comparison.png")
println("Saved: beam_deflection_stiffness_comparison.png")
display(p)
