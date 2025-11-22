include("plot_functions.jl")
include("mesh_geometry.jl")
using LinearAlgebra


######################
# STEP 1 create mesh # 
######################
n = 10 
L = 10 
H = 5 
sx = 1 
sy = 1 

points = generate_2Dmesh(L,H,n,sx,sy)

# plot_mesh(points)

# properties
K = 15        # Low stiffness for fluid elasticity
C = 10        # Higher damping for viscosity
total_mass = 2.2   # Slightly higher density

mass_nodes = total_mass/n^2 


# calculate natural length
natural_length_x = norm(points[1,2][1:2]-points[1,1][1:2])
natural_length_y = norm(points[2,1][1:2]-points[1,1][1:2])

function calculate_nodes!(i, j, dt)
    global natural_length_x, K, C, mass_nodes, natural_length_y
    
    # current node
    p = points[i, j]

    pos_p = p[1:2]
    vel_p = p[3:4]
    accn_p = p[5:6]
    
    # Initialize deltas and velocities with defaults for boundary handling
    del_east = 0
    del_west = 0
    del_north = 0
    del_south = 0
    vel_east = vel_p
    vel_west = vel_p
    vel_north = vel_p
    vel_south = vel_p

    # Get dimensions to check bounds
    rows, cols = size(points)
    
    # Check and calculate for adjacent nodes (only if in bounds)
    if i + 1 <= rows
       p_east = points[i+1, j][1:2]
       del_east = norm(p_east - pos_p) - natural_length_y
       vel_east = p[3:4]
    end
    
    if i - 1 >= 1
        p_west = points[i-1, j][1:2]
        del_west = norm(p_west - pos_p) - natural_length_y
        vel_west = p[3:4]
    end
    
    if j + 1 <= cols
        p_north = points[i, j+1][1:2]
        del_north = norm(p_north - pos_p) - natural_length_x
        vel_north = p[3:4]
    end
    
    if j - 1 >= 1
        p_south = points[i, j-1][1:2]
        del_south = norm(p_south - pos_p) - natural_length_x
        vel_south = p[3:4]
    end
    
    # Calculate new acceleration in x direction
    new_acc_x = (K * del_east - K*del_west - C* (vel_east[1] - vel_p[1]) - 
                 C * (vel_p[1] -  vel_west[1])) / mass_nodes
    
    # Calculate new acceleration in y direction
    new_acc_y = (K * del_north - K * del_south - C* (vel_north[2]
               - vel_p[2]) - C* (vel_p[2] - vel_south[2])) / mass_nodes
    
    # Update velocity and position
    new_vel_x = vel_p[1] + new_acc_x * dt
    new_pos_x = pos_p[1] + new_vel_x * dt
    
    new_vel_y = vel_p[2] + new_acc_y * dt
    new_pos_y = pos_p[2] + new_vel_y * dt
    
    # Update node state
    points[i, j] = [new_pos_x, new_pos_y, new_vel_x, new_vel_y, 
                    new_acc_x, new_acc_y]
end

function apply_fixed_left_boundary!()
    """
    Apply fixed boundary condition to the left edge (i=1).
    All nodes on the left edge (i=1) are fixed:
    - Position remains constant
    - Velocity set to zero
    - Acceleration set to zero
    """
    rows, cols = size(points)
    
    for j in 1:cols
        # Store the initial position
        x_fixed = points[1, j][1]
        y_fixed = points[1, j][2]
        
        # Set position, velocity, and acceleration
        points[1, j] = [x_fixed, y_fixed, 0.0, 0.0, 0.0, 0.0]
    end
end

# # Set initial conditions first
# points[end,1][1] = points[end,1][1] + 1
# points[end,1][2] = points[end,1][1] + 1
# points[end,1][3] = points[end,1][1] + 1

# Apply boundary condition to freeze left edge
apply_fixed_left_boundary!()

dt = 0.1
for i = 2:n, j = 1:n
    calculate_nodes!(i, j, dt)
end
