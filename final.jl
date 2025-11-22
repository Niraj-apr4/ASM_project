using LinearAlgebra
using Plots

# GEOMETRY
n_nodes = 11
nodes = [[Float64(i-1) * 0.1, 0.0] for i in 1:n_nodes]
edges = [[i, i+1] for i in 1:n_nodes-1]

# Material properties
k_spring = 1e4
damping = 0.9
B_field = [2.0, 1.0]
mag_moment = 20.0
fixed_nodes = [1]

# Simulation parameters
dt = 0.001
n_steps = 1000

# Initialize
pos = copy(nodes)
vel = [zeros(2) for _ in 1:n_nodes]
force = [zeros(2) for _ in 1:n_nodes]
mass = ones(n_nodes)

# Store position history
pos_history = [copy(pos)]

# Compute forces
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

# Run simulation
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

# Function to plot at discrete times
function plot_beam_at_times(pos_history, nodes, time_steps::Vector{Int}; B_field=[2.0, 1.0])
    plots_list = []
    
    for target_step in time_steps
        pos = pos_history[target_step]
        
        p = plot(xlim=(-0.1, 1.2), ylim=(-0.3, 0.5), 
                 aspect_ratio=:equal, legend=false,
                 title="Beam Deflection (Step $target_step)")
        
        # Original beam
        x_orig = [n[1] for n in nodes]
        y_orig = [n[2] for n in nodes]
        plot!(p, x_orig, y_orig, color=:gray, linestyle=:dash, linewidth=2)
        
        # Deformed beam
        x_def = [pt[1] for pt in pos]
        y_def = [pt[2] for pt in pos]
        plot!(p, x_def, y_def, color=:blue, linewidth=3)
        scatter!(p, x_def, y_def, color=:red, markersize=4)
        
        # Fixed node
        scatter!(p, [pos[1][1]], [pos[1][2]], 
                color=:green, markersize=8, markershape=:square)
        
        # Magnetic force arrow
        quiver!(p, [pos[end][1]], [pos[end][2]], 
                quiver=([B_field[1]*0.05], [B_field[2]*0.05]), 
                color=:purple, linewidth=2)
        
        push!(plots_list, p)
    end
    
    return plots_list
end

# Generate 3 plots at discrete times
plots = plot_beam_at_times(pos_history, nodes, [250, 500, 1000]; B_field=B_field)
combined = plot(plots[1], plots[2], plots[3], layout=(1,3), size=(1400, 400))
savefig(combined, "beam_deflection_3plots.png")
println("Saved: beam_deflection_3plots.png")
