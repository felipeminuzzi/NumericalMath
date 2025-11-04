import numpy as np
import matplotlib.pyplot as plt
import os

def plot_jato3D_results():
    """
    Reads the output .dat files from the jato3D Fortran program
    and generates contour and line plots.

    Assumes the .dat files are in the same directory as this script.
    
    The Fortran code writes 2D data files  and 1D centerline
    data files. This function plots both.
    """
    
    # --- 1. Define files for Contour Plots ---
    # Based on the Fortran write statements [cite: 36, 40]
    contour_files = {
        'u.dat': 'U Velocity (ua)',
        'v.dat': 'V Velocity (va)',
        'w.dat': 'W Velocity (wa)',
        'p.dat': 'Pressure (pa)',
        'iso.dat': 'Total Velocity Magnitude',
        'mixture.dat': 'Mixture Fraction (mista)',
        'viscosidade.dat': 'Effective Viscosity (mief)'
    }

    print("--- Generating Contour Plots ---")
    for filename, title in contour_files.items():
        if not os.path.exists(filename):
            print(f"Warning: File '{filename}' not found. Skipping plot.")
            continue
        
        try:
            # Load data. Skip the header row (e.g., "ni, nj, 1") 
            data = np.loadtxt(filename, skiprows=1)
            
            # Extract columns (x, y, variable) 
            x = data[:, 0]
            y = data[:, 1]
            var = data[:, 2]
            
            # Get unique coordinates to define the grid
            x_unique = np.unique(x)
            y_unique = np.unique(y)
            
            # Get grid dimensions
            ni = len(x_unique)
            nj = len(y_unique)
            
            # Create the 2D meshgrid
            X, Y = np.meshgrid(x_unique, y_unique)
            
            # Reshape the variable data to the (nj, ni) grid
            # The Fortran code loops j (rows) then i (cols),
            # so a standard (nj, ni) reshape works.
            VAR = var.reshape(nj, ni)
            
            # --- Create the Plot ---
            plt.figure(figsize=(12, 6))
            # Use 50 levels for a smooth contourf
            cp = plt.contourf(X, Y, VAR, levels=50, cmap='jet')
            plt.colorbar(cp, label=title)
            
            plt.title(f'Contour Plot of {title} (Z-Middle Plane)')
            plt.xlabel('X Coordinate')
            plt.ylabel('Y Coordinate')
            #plt.axis('equal') # Important for correct aspect ratio
            
            # Save the plot as a PNG
            plot_filename = f"{filename.split('.')[0]}_contour.png"
            plt.savefig(plot_filename)
            print(f"Saved: {plot_filename}")
            plt.show()

        except Exception as e:
            print(f"Error processing {filename}: {e}")

    # --- 2. Define files for Line Plots (Centerline) ---
    # Based on the Fortran write statements [cite: 36, 41]
    line_files = {
        'velcenter.dat': 'Centerline U Velocity',
        'mixcenter.dat': 'Centerline Mixture Fraction'
    }

    print("\n--- Generating Centerline Line Plots ---")
    for filename, title in line_files.items():
        if not os.path.exists(filename):
            print(f"Warning: File '{filename}' not found. Skipping plot.")
            continue

        try:
            # Load data. Skip the header row (e.g., "ni, 1") 
            data = np.loadtxt(filename, skiprows=1)
            
            # Extract columns (x, variable) 
            x = data[:, 0]
            var = data[:, 1]
            
            # --- Create the Plot ---
            plt.figure(figsize=(10, 6))
            plt.plot(x, var, marker='.', linestyle='-', markersize=4)
            
            plt.title(title)
            plt.xlabel('X Coordinate (Centerline)')
            plt.ylabel(title.split(' ')[-1]) # Use last word for label
            plt.grid(True, linestyle='--')
            
            # Save the plot as a PNG
            plot_filename = f"{filename.split('.')[0]}_line.png"
            plt.savefig(plot_filename)
            print(f"Saved: {plot_filename}")
            plt.show()

        except Exception as e:
            print(f"Error processing {filename}: {e}")

if __name__ == "__main__":
    # To run this script:
    # 1. Make sure you have numpy and matplotlib installed:
    #    pip install numpy matplotlib
    # 2. Run the Fortran code 'jato3D.f90' to generate the .dat files.
    # 3. Run this Python script in the same directory.
    plot_jato3D_results()