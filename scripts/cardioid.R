
# Define a range of values for t
t <- seq(0, 2 * pi, length.out = 1000)

# Calculate the x and y coordinates for the cardioid
x <- (1 - sin(t)) * cos(t)
y <- (1 - sin(t)) * sin(t)

x_A <- -0.6
y_A <- -1.0
theta_A <- atan2(y_A,x_A)
theta_search <- c(0,2*pi)

# plot.new()
# Create a plot of the cardioid
plot(x,y, type = "l", col = "blue", xlab = "X", ylab = "Y", main = "Cardioid")
points(x_A,y_A, col = "red", pch = 19)
segments(0, 0, x_A, y_A, col = "red", lwd = 2)

min_distance <- Inf
min_x <- Inf
min_y <- Inf
for (t2 in seq(theta_search[1], theta_search[2], by = 0.001)) {
  x_coord <- (1 - sin(t2)) * cos(t2)
  y_coord <- (1 - sin(t2)) * sin(t2)
  
  distance <- sqrt((x_coord - x_A)^2 + (y_coord - y_A)^2)
  if (distance < min_distance) {
    min_distance <- distance
    min_x <- x_coord
    min_y <- y_coord
  }
}
segments(min_x, min_y, x_A, y_A, col = "green", lwd = 2)
# Add a grid for reference (optional)
grid()