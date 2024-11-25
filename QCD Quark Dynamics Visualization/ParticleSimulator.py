# Re-import necessary libraries after reset
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd   
from matplotlib.gridspec import GridSpec
from scipy.optimize import fsolve

# Quark colors
quark_colors = ['red', 'blue', 'green']
 
# Constants
alpha_s = 0.2  # Approximate value of the strong coupling constant
sigma = 0.9  # String tension (GeV/fm, approximate value)
quark_positions_3d = np.array([
    [0, 0, 0],      # Quark 1
    [0.6, 0, 0],      # Quark 2
    [-0.1, 0.6, 0.6]  # Quark 3
])

# Function to calculate the potential energy and force between two quarks
def qcd_potential_force(r, alpha_s, sigma):
    if r == 0:
        return 0, 0  # Avoid division by zero
    potential = -4/3 * alpha_s / r + sigma * r  # QCD potential
    force_magnitude = -(-4/3 * alpha_s / r**2 + sigma)  # Derivative of potential
    return potential, force_magnitude 
 

# Initialize lists to keep track of scatter plot objects
quark_potential_scatters = []
quark_force_scatters = []


# Update function for animation with 2D plot updates
def update(frame):
    # Update 3D plot (same as before)
    current_positions = positions_over_time[frame]
    quark_scatter._offsets3d = (
        current_positions[:, 0],
        current_positions[:, 1],
        current_positions[:, 2]
    )
    quark_scatter.set_color(quark_colors)
    
    # Update lines and force labels in 3D plot
    for i, line in enumerate(lines):
        p1 = current_positions[i]
        p2 = current_positions[(i + 1) % len(current_positions)]
        line.set_data([p1[0], p2[0]], [p1[1], p2[1]])
        line.set_3d_properties([p1[2], p2[2]])
        
        # Calculate force between connected quarks
        r_vec = p2 - p1
        r_mag = np.linalg.norm(r_vec)
        _, force_mag = qcd_potential_force(r_mag, alpha_s, sigma)
        #line.set_color(force_color(force_mag))  # Change line color based on force magnitude
        
        # Update force label
        midpoint = (p1 + p2) / 2
        force_labels[i].set_position((midpoint[0], midpoint[1]))
        force_labels[i].set_3d_properties(midpoint[2])
        force_labels[i].set_horizontalalignment('center')  # Ensure horizontal alignment
        force_labels[i].set_text(f"{force_mag:.2f} GeV/fm")
    
    # Update 2D plots for quark potentials and forces
    distances = [
        np.linalg.norm(current_positions[i] - current_positions[j])
        for i in range(len(current_positions))
        for j in range(i + 1, len(current_positions))
    ]
    
    # Remove previous scatter plots
    for scatter in quark_potential_scatters:
        scatter.remove()
    for scatter in quark_force_scatters:
        scatter.remove()
    
    # Clear the lists
    quark_potential_scatters.clear()
    quark_force_scatters.clear()
    
    quark_potentials = []
    quark_forces = []
    for r in distances:
        potential, force_mag = qcd_potential_force(r, alpha_s, sigma)
        quark_potentials.append(potential)
        quark_forces.append(force_mag)
    
    # Add scatter points for updated quark potentials and forces
    quark_potential_scatters.append(ax2.scatter(
        distances,
        quark_potentials,
        color='orange',
        s=100,
        zorder=5,
        label='Quark Potentials' if frame == 0 else ""
    ))
    quark_force_scatters.append(ax2_f.scatter(
        distances,
        quark_forces,
        color='purple',
        s=100,
        zorder=5,
        label='Quark Forces' if frame == 0 else ""
    )) 
    
    # Update legends dynamically
    ax2.legend(loc='upper left')
    ax2_f.legend(loc='upper right')

    return [quark_scatter, *lines, *force_labels, *quark_potential_scatters, *quark_force_scatters]


# Define a function for the force derived from the QCD potential
def force_function(r, alpha_s=0.2, sigma=0.9):
    _, force = qcd_potential_force(r, alpha_s, sigma)
    return force
 


from matplotlib.animation import FuncAnimation

# Time parameters for simulation
time_steps = 1000  # Number of time steps
delta_t = 0.1  # Time step duration (arbitrary units)

# Initialize positions and velocities
positions = quark_positions_3d.copy()
velocities = np.zeros_like(positions)

# Store positions over time for visualization
positions_over_time = [positions.copy()]

# Simulate quark motion
for _ in range(time_steps):
    forces = np.zeros_like(positions)
    for i in range(len(positions)):
        for j in range(len(positions)):
            if i != j:
                r_vec = positions[j] - positions[i]
                r_mag = np.linalg.norm(r_vec)
                _, force_mag = qcd_potential_force(r_mag, alpha_s, sigma)
                force_dir = r_vec / r_mag if r_mag != 0 else 0
                forces[i] -= force_mag * force_dir
                
    # Update velocities and positions
    velocities += forces * delta_t  # Simplified dynamics
    positions += velocities * delta_t
    positions_over_time.append(positions.copy())

# Convert list to numpy array for easier handling
positions_over_time = np.array(positions_over_time)

 
 
# QCD Potential and Force vs. Distance Graph
  
# Define a range of distances (fm) for plotting
r_values = np.linspace(0.01, 2.0, 500)  # Avoid zero to prevent division by zero
potentials = []
forces = []

# Calculate QCD potential and force for each distance
for r in r_values:
    potential, force_mag = qcd_potential_force(r, alpha_s, sigma)
    potentials.append(potential)
    forces.append(force_mag)

# Solve for the equilibrium point where the force is zero
initial_guess = 1.0  # Start near typical quark separation (fm)
equilibrium_position = fsolve(force_function, initial_guess)[0]



# figure
fig = plt.figure(figsize=(20, 8))
gs = GridSpec(1, 2, width_ratios=[1, 1])


# Create 3D animation
ax1 = fig.add_subplot(gs[0], projection='3d')
ax1.set_title('Quark Dynamics in a Proton (3D Animation)')

# Initialize plot elements
quark_scatter = ax1.scatter([], [], [], color='red', s=100, label='Quarks')
lines = [ax1.plot([], [], [], 'k--', alpha=0.7)[0] for _ in range(3)]
force_labels = [ax1.text(0, 0, 0, '', color='blue') for _ in range(3)]

# Set up plot limits and labels
ax1.set_xlim(np.min(positions_over_time), np.max(positions_over_time))
ax1.set_ylim(np.min(positions_over_time), np.max(positions_over_time))
ax1.set_zlim(np.min(positions_over_time), np.max(positions_over_time))
ax1.set_xlabel('x (fm)')
ax1.set_ylabel('y (fm)')
ax1.set_zlabel('z (fm)')
ax1.legend() 


# Graph 2: QCD Potential and Force
ax2 = fig.add_subplot(gs[1])
ax2.set_title('QCD Potential and Force vs. Distance (Equilibrium Highlighted)', fontsize=14)

# Plot potential energy
ax2.plot(r_values, potentials, 'r-', label='Potential (GeV)')
ax2.set_xlabel('Distance (fm)', fontsize=12)
ax2.set_ylabel('Potential (GeV)', color='r', fontsize=12)
ax2.tick_params(axis='y', labelcolor='r')
ax2.axhline(0, color='gray', linestyle='--', linewidth=0.7, label='Potential = 0')
ax2.scatter(
    [equilibrium_position],
    [qcd_potential_force(equilibrium_position, alpha_s, sigma)[0]],
    color='red',
    label='Equilibrium Point',
    zorder=5
)

# Plot force on secondary axis
ax2_f = ax2.twinx()
ax2_f.plot(r_values, forces, 'b-', label='Force (GeV/fm)')
ax2_f.set_ylabel('Force (GeV/fm)', color='b', fontsize=12)
ax2_f.tick_params(axis='y', labelcolor='b')
ax2_f.axhline(0, color='gray', linestyle='--', linewidth=0.7, label='Force = 0')
ax2_f.scatter(
    [equilibrium_position],
    [qcd_potential_force(equilibrium_position, alpha_s, sigma)[0]],
    color='blue',
    label='Equilibrium Point',
    zorder=5
)


fig.tight_layout()
  
# Add combined legend
ax2.legend(loc='upper left')
ax2_f.legend(loc='upper right')

# Create animation
ani = FuncAnimation(fig, update, frames=time_steps, interval=50, blit=False)

# Display the side-by-side plots
plt.show()

# Display equilibrium position
equilibrium_position

