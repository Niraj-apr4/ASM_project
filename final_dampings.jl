using LinearAlgebra
using Plots

# GEOMETRY
n_nodes = 11
nodes = [[Float64(i-1) * 0.1, 0.0] for i in 1:n_nodes]
edges = [[i, i+1] for i in 1:n_nodes-1]

# Material properties
k_spring = 1e4
mag_moment = 50.0  # Fixed magnetic moment
B_field = [1.0, 1.0]
damping_coefficients = [0.5, 0.7, 1.0]
fixed_nodes = [1]

# Simulation parameters
dt = 0.01
n_steps = 1000

# Function to run simulation for a given damping coefficient
function run_simulation(damping, n_steps, dt)
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

# Run simulations for all damping coefficients
all_pos_histories = []
for damping_coeff in damping_coefficients
    println("Running simulation for damping = $damping_coeff")
    pos_hist = run_simulation(damping_coeff, n_steps, dt)
    push!(all_pos_histories, pos_hist)
end

# Create subplots for each damping coefficient
final_step = n_steps + 1
subplots_list = []

for (idx, damping_coeff) in enumerate(damping_coefficients)
    pos_hist = all_pos_histories[idx]
    final_pos = pos_hist[final_step]
    
    p = plot(xlim=(-0.1, 1.2), ylim=(-1.2, 0.5), 
             aspect_ratio=:equal, legend=false,
             title="Damping Coefficient: c = $damping_coeff",
             xlabel="x (m)", ylabel="y (m)",
             size=(400, 400))
    
    # Plot original beam
    x_orig = [n[1] for n in nodes]
    y_orig = [n[2] for n in nodes]
    plot!(p, x_orig, y_orig, color=:gray, linestyle=:dash, linewidth=2, alpha=0.5)
    
    # Plot deformed beam
    x_def = [pt[1] for pt in final_pos]
    y_def = [pt[2] for pt in final_pos]
    
    plot!(p, x_def, y_def, color=:blue, linewidth=3, alpha=0.8)
    scatter!(p, x_def, y_def, color=:red, markersize=4, alpha=0.7)
    
    # Highlight fixed node
    scatter!(p, [final_pos[1][1]], [final_pos[1][2]], 
             color=:green, markersize=8, markershape=:square)
    
    # Magnetic force arrow
    quiver!(p, [final_pos[end][1]], [final_pos[end][2]], 
            quiver=([B_field[1]*0.05], [B_field[2]*0.05]), 
            color=:purple, linewidth=2)
    
    push!(subplots_list, p)
end

# Create combined figure with 3 subplots
combined = plot(subplots_list[1], subplots_list[2], subplots_list[3], 
                layout=(1,3), size=(1400, 400),
                plot_title="Beam Deflection for Different Damping Coefficients (μ = $(mag_moment))")

savefig(combined, "beam_deflection_damping_comparison.png")
println("Saved: beam_deflection_damping_comparison.png")
display(combined)
