#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 15:01:29 2024

@author: micahpratt
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import numpy as np
import matplotlib.pyplot as plt
import os
import cv2
import imageio
import math

# Load the data from the file
data_path = '/Users/micahpratt/Documents/LAB/Landslide/DEM/multi_trajectory.txt'
if not os.path.exists(data_path):
    raise FileNotFoundError(f"Data file not found at: {data_path}")

data = np.loadtxt(data_path, comments='#')

# Extract time, x, and y positions
npoints = 64


total_frames = len(data[:, 1])
start_Frame = 0
end_Frame = total_frames # round(total_frames * 3/8)

omega = data[start_Frame:end_Frame, 7]
theta = data[start_Frame:end_Frame, 6]
radius = data[start_Frame:end_Frame, 1] 
time = data[start_Frame:end_Frame, 0]
x_positions = data[start_Frame:end_Frame, 2]
y_positions = data[start_Frame:end_Frame , 3]



# Create a folder to save frames
frames_dir = "/Users/micahpratt/Documents/LAB/DEM/frames"
if not os.path.exists(frames_dir):
    os.makedirs(frames_dir)

# Plot and save each frame

for i in range(int(len(time) / npoints)):
    fig = plt.figure(figsize=(6.4, 4.8))
    ax = fig.add_subplot()

    # Define start and end rows for the current frame
    start_row = i * npoints
    end_row = start_row + npoints

    x_positions_step = x_positions[start_row:end_row]
    y_positions_step = y_positions[start_row:end_row]
    radius_step = radius[start_row:end_row]
    omega_step = omega[start_row:end_row]
    theta_step = theta[start_row:end_row]

    # Loop through each ball in this frame
    for k in range(npoints):
        color = 'k'
        # Determine circle color based on omega
        if omega_step[k] < 1e-6 and omega_step[k] > 0:
            color = 'k'  # Black for near-zero omega
        elif omega_step[k] < 0:
            color = 'b'  # Blue for negative omega
        elif omega_step[k] > 0:
            color = 'r'  # Red for positive omega

        # Create the circle
        circle = plt.Circle((x_positions_step[k], y_positions_step[k]), radius_step[k], color=color)
        ax.add_patch(circle)

        # Add white rotation line
        xline = np.array([
            x_positions_step[k],
            x_positions_step[k] + radius_step[k] * math.cos(theta_step[k])
        ])
        yline = np.array([
            y_positions_step[k],
            y_positions_step[k] + radius_step[k] * math.sin(theta_step[k])
        ])
        plt.plot(xline, yline, "w-", lw=2.0)  # Draw the rotation line

    # Set plot limits and labels
    plt.xlim(-2, 2)
    plt.ylim(0, 2)
    plt.xlabel("X Position (m)")
    plt.ylabel("Y Position (m)")
    plt.title(f"Time: {time[i * npoints]:.3f}s")
    ax.set_aspect('equal')

    # Save the frame as an image
    frame_filename = os.path.join(frames_dir, f"frame_{i:04d}.png")
    plt.savefig(frame_filename, dpi=300)
    plt.close()
    print(f"Saved frame {i}: {frame_filename}")

# Parameters for the video
frame_rate = 30
output_gif_path = '/Users/micahpratt/Documents/LAB/DEM/bouncing_ball.gif'
images = []

# Read and append frames for GIF creation
for i in range(int(len(time) / npoints)):
    frame_filename = os.path.join(frames_dir, f"frame_{i:04d}.png")
    if os.path.exists(frame_filename):
        images.append(imageio.imread(frame_filename))
    else:
        print(f"Missing frame: {frame_filename}")

# Create GIF
kargs = {'duration': 1 / frame_rate}
imageio.mimsave(output_gif_path, images, 'GIF', **kargs)
print(f"GIF saved at: {output_gif_path}")
