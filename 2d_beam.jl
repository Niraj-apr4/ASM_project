using LinearAlgebra
using Plots

# GEOMETRY
n_x = 11
n_y = 11
n_nodes = n_x * n_y

# Create nodes in 2D grid
nodes = []
for j in 1:n_y
    for i in 1:n_x
        x = Float64(i-1) * 0.1
        y = Float64(j-1) * 0.1
        push!(nodes, [x, y])
    end
end

# Function to get node index from grid position
function get_index(i, j)
    return (j-1) * n_x + i
end

# Create edges (spring connections between neighbors)
edges = []
for j in 1:n_y
    for i in 1:n_x
        # Right neighbor
        if i < n_x
            push!(edges, [get_index(i, j), get_index(i+1, j)])
        end
        # Bottom neighbor
        if j < n_y
            push!(edges, [get_index(i, j), get_index(i, j+1)])
        end
    end
end

# Material properties
k_spring = 1e4
damping = 0.9
B_field = [1.0, 1.0]
mag_moment = 50.0
fixed_nodes = [get_index(i, j) for i in 1:1 for j in 1:n_y]  # Left edge (x=0)

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
    
    # Spring forces
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
    
    # Magnetic force on free boundary nodes (right edge, x=1.0)
    for j in 1:n_y
        idx = get_index(n_x, j)
        force[idx] += mag_moment * B_field
    end
    
    # Gravitational force
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

# Plotting function
function plot_2d_beam(pos, nodes, title_str="2D Beam Deflection")
    p = plot(xlim=(-0.2, 1.3), ylim=(-0.3, 1.3), 
             aspect_ratio=:equal, legend=false,
             title=title_str, xlabel="x (m)", ylabel="y (m)",
             size=(700, 700))
    
    # Plot original grid (faint)
    for edge in edges
        i, j = edge
        x_orig = [nodes[i][1], nodes[j][1]]
        y_orig = [nodes[i][2], nodes[j][2]]
        plot!(p, x_orig, y_orig, color=:gray, linestyle=:dash, linewidth=0.5, alpha=0.3)
    end
    
    # Plot deformed grid
    for edge in edges
        i, j = edge
        x_def = [pos[i][1], pos[j][1]]
        y_def = [pos[i][2], pos[j][2]]
        plot!(p, x_def, y_def, color=:blue, linewidth=1, alpha=0.6)
    end
    
    # Plot nodes
    x_nodes = [pt[1] for pt in pos]
    y_nodes = [pt[2] for pt in pos]
    scatter!(p, x_nodes, y_nodes, color=:red, markersize=3, alpha=0.5)
    
    # Highlight fixed nodes (left edge)
    for j in 1:n_y
        idx = get_index(1, j)
        scatter!(p, [pos[idx][1]], [pos[idx][2]], 
                color=:green, markersize=5, markershape=:square)
    end
    
    # Highlight loaded nodes (right edge)
    for j in 1:n_y
        idx = get_index(n_x, j)
        scatter!(p, [pos[idx][1]], [pos[idx][2]], 
                color=:purple, markersize=5)
    end
    
    return p
end

# Generate plots at different timesteps
final_step = n_steps + 1
p1 = plot_2d_beam(pos_history[250], nodes, "2D Beam Deflection (Step 250)")
p2 = plot_2d_beam(pos_history[500], nodes, "2D Beam Deflection (Step 500)")
p3 = plot_2d_beam(pos_history[final_step], nodes, "2D Beam Deflection (Step 1000 - Equilibrium)")

# Combine plots
combined = plot(p1, p2, p3, layout=(1,3), size=(1800, 600))
savefig(combined, "beam_2d_deflection.png")
println("Saved: beam_2d_deflection.png")
display(combined)
