"""
contour_plot.jl
objectives : for plotting contour plot with different
            mesh geometry and precision
"""

using Plots

function plot_mesh(nodes)
    points = vec(nodes)
    x = [p[1] for p in points]
    y = [p[2] for p in points]
    
    # Get min and max values
    x_min, x_max = minimum(x), maximum(x)
    y_min, y_max = minimum(y), maximum(y)
    
    s = scatter(x, y, 
            aspect_ratio=:equal,
            markersize=1,
            xlabel="X",
            ylabel="Y",
            legend=false,
            title="Computational Nodes",
            xticks=[x_min, x_max],
            yticks=[y_min, y_max],
            fontfamily="Computer Modern")
    display(s)
end

"""
TODO  

"""
function plot_temp(data)

    # Get unique coordinates
    x_vals = sort(unique([p[1] for p in data]))
    y_vals = sort(unique([p[2] for p in data]))
    
    # Create temperature matrix
    temp_matrix = zeros(length(y_vals), length(x_vals))
    for point in data
        i = findfirst(y -> y ≈ point[2], y_vals)
        j = findfirst(x -> x ≈ point[1], x_vals)
        temp_matrix[i, j] = point[3]
    end
    
    # contour plot
    c = contourf(x_vals, y_vals, temp_matrix,
            title="Temperature Distribution",
            xlabel="X Coordinate (m)", 
            ylabel="Y Coordinate (m)",
            colorbar_title="Temperature (°C)",
                 clabels = true,
                 linecolor=:white,
            aspect_ratio=:equal,
            levels=25,
            color=:turbo,
            size=(700, 600),
            dpi=300,
            linewidth=0,
            fontfamily="Computer Modern")
    display(c)
end

function plot_velocity(nodes)
    
    # Extract x, y, u, v from nodes
    n_i, n_j = size(nodes)
    xp = [nodes[i,j][1] for i in 1:n_i, j in 1:n_j]
    yp = [nodes[i,j][2] for i in 1:n_i, j in 1:n_j]
    u = [nodes[i,j][4] for i in 1:n_i, j in 1:n_j]
    v = [nodes[i,j][5] for i in 1:n_i, j in 1:n_j]
    
    # Plot velocity field
    vec = 5
    q = quiver(xp, yp, quiver=(u', v', vec), 
           aspect_ratio=:equal, 
           xlabel="x", ylabel="y",
           legend=false)
    display(q)
end

# Plot with log scaling on x-axis (error tolerance)
function plot_error(error_trial, counter_array)
    error_vec = vec(error_trial)
    counter_vec = vec(counter_array)
    plt = plot(
        error_vec, counter_vec;
        xscale = :log10,
        xlabel = "Error tolerance",
        ylabel = "Number of iterations",
        # TODO implement a variable way to change the title of
        # plot
        title = " TDMA Convergence Test",
        legend = false,
        lw = 2,
        marker = :circle,
        markersize = 6,
        markercolor = :blue,
        linecolor = :blue,
        grid = :on,
        framestyle = :box,
        background_color = :white,
        guidefont = font(12, "Computer Modern"),
        tickfont = font(10, "Computer Modern"),
        titlefont = font(14, "Computer Modern")
    )

    display(plt)
    return plt
end

function plot_mesh_grid(points, title_str="Mesh Grid")
    rows, cols = size(points)
    
    # Extract positions and velocities
    positions = zeros(2, rows, cols)
    velocities = zeros(2, rows, cols)
    
    for i in 1:rows, j in 1:cols
        positions[1, i, j] = points[i, j][1]
        positions[2, i, j] = points[i, j][2]
        velocities[1, i, j] = points[i, j][3]
        velocities[2, i, j] = points[i, j][4]
    end
    
    # Create plot
    p = plot(legend=false, aspect_ratio=:equal, title=title_str)
    
    # Plot grid lines (horizontal)
    for i in 1:rows
        x_vals = positions[1, i, :]
        y_vals = positions[2, i, :]
        plot!(p, x_vals, y_vals, linewidth=2, color=:blue, alpha=0.7)
    end
    
    # Plot grid lines (vertical)
    for j in 1:cols
        x_vals = positions[1, :, j]
        y_vals = positions[2, :, j]
        plot!(p, x_vals, y_vals, linewidth=2, color=:blue, alpha=0.7)
    end
    
    # Plot nodes as scatter points
    x_nodes = vec(positions[1, :, :])
    y_nodes = vec(positions[2, :, :])
    scatter!(p, x_nodes, y_nodes, markersize=6, color=:red, alpha=0.8)
    
    xlabel!("X")
    ylabel!("Y")
    return p
end

function create_beam_animation(num_steps, frame_interval, filename="magnetoelastic_beam_animation.mp4")
    """
    Create animation of magnetoelastic beam simulation
    
    Parameters:
    - num_steps: total simulation steps
    - frame_interval: save frame every N steps
    - filename: output MP4 filename
    """
    
    anim = @animate for step in 1:num_steps
        # Update nodes
        for i = 2:n, j = 1:n
            calculate_nodes!(i, j, dt)
        end
        apply_fixed_left_boundary!()
        
        # Save frame every frame_interval steps
        if mod(step, frame_interval) == 0
            rows, cols = size(points)
            positions = zeros(2, rows, cols)
            
            for i in 1:rows, j in 1:cols
                positions[1, i, j] = points[i, j][1]
                positions[2, i, j] = points[i, j][2]
            end
            
            p = plot(legend=false, aspect_ratio=:equal, 
                     title="Magnetoelastic Fluid Beam - Step: $step",
                     xlim=(-0.5, L+0.5), ylim=(-0.5, H+0.5), size=(1400, 900))
            
            # Plot horizontal grid lines
            for i in 1:rows
                plot!(p, positions[1, i, :], positions[2, i, :], 
                      linewidth=2, color=:blue, alpha=0.7)
            end
            
            # Plot vertical grid lines
            for j in 1:cols
                plot!(p, positions[1, :, j], positions[2, :, j], 
                      linewidth=2, color=:blue, alpha=0.7)
            end
            
            # Plot nodes
            scatter!(p, vec(positions[1, :, :]), vec(positions[2, :, :]), 
                     markersize=6, color=:red, alpha=0.8)
            
            # Highlight fixed boundary
            scatter!(p, positions[1, 1, :], positions[2, 1, :], 
                     markersize=8, color=:green, alpha=0.9)
            
            xlabel!("X")
            ylabel!("Y")
        end
    end
    
    # Save animation
    mp4(anim, filename, fps=20)
    println("Animation saved as '$filename'")
end
