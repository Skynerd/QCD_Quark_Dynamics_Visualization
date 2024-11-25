# QCD Quark Dynamics Visualization

## Overview

This project simulates and visualizes the dynamics of quarks in a proton, governed by Quantum Chromodynamics (QCD) principles. It incorporates:

- **Quark dynamics simulation**: Models the motion of quarks under QCD potential and force.
- **QCD potential and force analysis**: Illustrates the interplay of forces and energy between quarks as a function of distance.
- **3D animations and graphs**: Provides interactive and static visualizations for intuitive understanding of QCD concepts.

The project highlights the equilibrium position where the net force between quarks vanishes, a critical concept in understanding proton stability.

---

## Features

- **3D Simulation**: Animates the positions of three quarks interacting under QCD forces.
- **Dynamic Force Visualization**: Displays time-dependent force magnitudes between quarks.
- **Potential and Force Plots**: Plots QCD potential and force as functions of distance, with equilibrium points marked.
- **Equilibrium Calculation**: Solves for the distance at which forces balance using numerical techniques.

---

## Requirements

To run this project, install the following Python packages:

- `numpy`
- `matplotlib`
- `scipy`
- `pandas`

---

## Code Highlights

### Core Functions

- **`qcd_potential_force(r, alpha_s, sigma)`**: Computes the QCD potential and force magnitude for a given quark distance.
- **`force_function(r)`**: Derives the force from the QCD potential, used to find equilibrium points.
- **`update(frame)`**: Updates the 3D animation and dynamic plots for each frame.

### Simulation

- Quark positions evolve over time, considering forces derived from the QCD potential.
- Numerical integration updates quark positions and velocities.

### Visualization

1. **3D Animation**: 
   - Quark positions, force vectors, and dynamic labels.
   - Line connections between quarks, updated in real time.

2. **QCD Potential and Force Plots**:
   - QCD potential (`V(r) = -4/3 * α_s / r + σ * r`).
   - Force derived from potential.
   - Highlights equilibrium points where force is zero.

---

## How to Use

1. Clone the repository and navigate to the project directory.
2. Run the script to visualize the QCD dynamics:
   ```bash
   python qcd_quark_dynamics.py
   ```
3. The following visualizations will appear:
   - **3D animation** of quark dynamics.
   - **Graphs** showing QCD potential and force vs. distance, including equilibrium points.

---

## Example Output

1. **3D Animation**:
   - Shows quark motion influenced by QCD forces.
   - Updates quark positions and visualizes forces dynamically.

2. **Graphs**:
   - **Potential vs. Distance**: Illustrates the energy landscape of quark interactions.
   - **Force vs. Distance**: Highlights where the net force vanishes (equilibrium).

---

## Key Physics Concepts

- **Strong Coupling Constant (αₛ)**: Determines the strength of the strong force.
- **String Tension (σ)**: Represents the confinement of quarks due to gluonic flux tubes.
- **QCD Potential**:
  \[
  V(r) = -\frac{4}{3} \frac{\alpha_s}{r} + \sigma r
  \]
- **Equilibrium**: The distance where attractive and repulsive forces between quarks balance.

---

## Contact

For questions or collaboration, feel free to reach out or submit issues in the repository.
