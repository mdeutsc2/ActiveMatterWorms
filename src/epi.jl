using Plots
using LinearAlgebra

# Parametric equations for the epitrochoid
function epitrochoid_points1(R, r, d, spacing)
    # Initialize arrays to store the points
    points = []

    # Function to calculate x, y coordinates for a given angle θ
    function epitrochoid_coords(θ)
        x = (R + r) * cos(θ) - d * cos((R + r) / r * θ)
        y = (R + r) * sin(θ) - d * sin((R + r) / r * θ)
        return [x, y]
    end

    # Arc length estimation function (using numerical integration)
    function arc_length(θ1, θ2)
        n = 10000  # Number of subdivisions for numerical integration
        deltaθ = (θ2 - θ1) / n
        length = 0.0
        for i in 1:n-1
            θ_left = θ1 + (i - 1) * deltaθ
            θ_right = θ1 + i * deltaθ
            p_left = epitrochoid_coords(θ_left)
            p_right = epitrochoid_coords(θ_right)
            # Use Euclidean distance between consecutive points
            length += norm(p_right - p_left)
        end
        return length
    end

    # We will start at θ = 0 and loop while we are within the desired arc length
    θ = 0.0
    total_length = 0.0
    # Store the first point
    push!(points, epitrochoid_coords(θ))
    println("Est: # pts: ",floor(arc_length(0.0, 2 * π)/spacing))
    
    while total_length < arc_length(0.0, 2 * π)  # Loop through one full revolution
        # Estimate the next point by stepping along θ
        step_size = 0.001  # Small step to refine placement
        next_θ = θ + step_size
        
        # Compute the arc length increment
        current_length = arc_length(0.0, next_θ)
        
        # If the arc length difference is greater than or equal to the spacing, place a point
        if current_length - total_length >= spacing
            # Append point
            push!(points, epitrochoid_coords(next_θ))
            total_length = current_length  # Update the total length
            θ = next_θ  # Move to the next step
        else
            # Otherwise just increment θ by the small step
            θ = next_θ
        end
    end

    return points
end

# Parametric equations for the epitrochoid
function epitrochoid_points2(R, r, d, num_points)
    # Initialize array to store the points
    points = []

    # Function to calculate x, y coordinates for a given angle θ
    function epitrochoid_coords(θ)
        x = (R + r) * cos(θ) - d * cos((R + r) / r * θ)
        y = (R + r) * sin(θ) - d * sin((R + r) / r * θ)
        return [x, y]
    end

    # Arc length estimation function (using numerical integration)
    function arc_length(θ1, θ2)
        n = 10000  # Number of subdivisions for numerical integration
        deltaθ = (θ2 - θ1) / n
        length = 0.0
        for i in 1:n-1
            θ_left = θ1 + (i - 1) * deltaθ
            θ_right = θ1 + i * deltaθ
            p_left = epitrochoid_coords(θ_left)
            p_right = epitrochoid_coords(θ_right)
            # Use Euclidean distance between consecutive points
            length += norm(p_right - p_left)
        end
        return length
    end

    # Calculate the total arc length for one full revolution (0 to 2π)
    total_length = arc_length(0.0, 2 * π)

    # Calculate the spacing between points
    spacing = total_length / (num_points - 1)

    # Initialize the first point
    θ = 0.0
    total_length_travelled = 0.0
    push!(points, epitrochoid_coords(θ))

    # Loop through and place evenly spaced points
    for _ in 2:num_points
        while true
            # Increment angle by a small step
            θ += 0.001
            current_length = arc_length(0.0, θ)
            # Check if we have traveled enough to place the next point
            if current_length - total_length_travelled >= spacing
                push!(points, epitrochoid_coords(θ))
                total_length_travelled = current_length  # Update total length traveled
                break
            end
        end
    end

    return points
end

function epitrochoid_points3(R,r,d,num_points)
    # Initialize array to store the points
    points = []

    # Function to calculate x, y coordinates for a given angle θ
    function epitrochoid_coords(θ)
        x = (R + r) * cos(θ) - d * cos((R + r) / r * θ)
        y = (R + r) * sin(θ) - d * sin((R + r) / r * θ)
        return [x, y]
    end

    # Loop through and place evenly spaced points
    for itheta in range(start=0,stop=2*pi,length=num_points)
        push!(points, epitrochoid_coords(itheta))
    end
    return points
end

# Example usage:
k = 2 # number of cusps
r = 2.0  # Radius of the rolling circle
R = k*r # Radius of the fixed circle
println(r)
println(R)
d = 0.9999  # Distance of the point from the center of the rolling circle
num_points = 100  # Desired number of evenly spaced points

# Generate points along the perimeter of the epitrochoid
#points = epitrochoid_points2(R, r, d, num_points)
points = epitrochoid_points1(R,r,d,0.15)
#points = epitrochoid_points3(R, r, d, num_points)
println(length(points))
for i in 1:10
    println(points[i])
end
# Extract x and y coordinates from the points
x_vals = [p[1] for p in points]
y_vals = [p[2] for p in points]

# Plot the points
scatter(x_vals, y_vals, label="Epitrochoid", xlabel="x", ylabel="y", legend=:topright, aspect_ratio=:equal, show=true,markersize=0.25)
