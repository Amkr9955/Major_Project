import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path
import matplotlib as mpl

# Set research-quality plotting parameters
mpl.rcParams['figure.dpi'] = 300
mpl.rcParams['font.size'] = 12
mpl.rcParams['axes.labelsize'] = 14
mpl.rcParams['axes.titlesize'] = 16
mpl.rcParams['xtick.labelsize'] = 12
mpl.rcParams['ytick.labelsize'] = 12
mpl.rcParams['legend.fontsize'] = 12
mpl.rcParams['figure.autolayout'] = True

def plot_data(T, Y, y_label, filename, T_min=None, T_max=None, T_step=None, is_d_plot=False):
    # Create mask for temperature range
    mask = np.ones_like(T, dtype=bool)
    if T_min is not None:
        mask = mask & (T >= T_min)
    if T_max is not None:
        mask = mask & (T <= T_max)
    
    # Apply temperature range filter
    T_filtered = T[mask]
    Y_filtered = Y[mask]
    
    # Apply temperature step filter if specified
    if T_step is not None and len(T_filtered) > 0:
        # Find unique temperature values within the range
        unique_T = np.unique(T_filtered)
        # Select temperatures at the specified step
        step = int(T_step / np.min(np.diff(unique_T)))
        selected_T = unique_T[::step] if step > 0 else unique_T
        # Create mask for selected temperatures
        step_mask = np.isin(T_filtered, selected_T)
        T_filtered = T_filtered[step_mask]
        Y_filtered = Y_filtered[step_mask]
    
    # Sort data for proper visualization
    sorted_indices = np.argsort(T_filtered)
    T_sorted, Y_sorted = T_filtered[sorted_indices], Y_filtered[sorted_indices]

    # Customize plot parameters for d1, d2, d3 plots
    if is_d_plot:
        plt.figure(figsize=(10, 7))  # Larger figure size
        plt.plot(T_sorted, Y_sorted, color='black', linestyle="-", marker='o', 
                markersize=4, linewidth=1.5, label='Data')  # Bigger markers and thicker lines
        plt.xlabel('Temperature (T)', fontsize=14)
        plt.ylabel(y_label, fontsize=14)
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(fontsize=12)
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.title(f'{y_label} vs Temperature', fontsize=16)
    else:
        plt.figure(figsize=(7, 5))
        plt.plot(T_sorted, Y_sorted, color='black', linestyle="-", marker='o', 
                markersize=3, label='Data')
        plt.xlabel('Temperature (T)')
        plt.ylabel(y_label)
        plt.legend()
        plt.grid()
        plt.title(f'{y_label} vs Temperature')
    
    plt.tight_layout()  # Adjust layout to prevent label cutoff
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()

def create_combined_plots(all_data, output_dir):
    """Create combined plots for Cv, C1, C2, and C3"""
    # Define colors and markers for different files
    colors = plt.cm.viridis(np.linspace(0, 1, len(all_data)))
    markers = ['o', 's', '^', 'D', 'v', 'p', '*', 'h', '8', '<', '>']
    
    # Create combined Cv plot
    plt.figure(figsize=(15, 12))
    for i, (file_stem, data) in enumerate(all_data.items()):
        T = data["Temperature"].values
        Cv = data["Cv"].values
        plt.plot(T, Cv, linestyle='-', marker=markers[i % len(markers)], 
                markersize=5, linewidth=1.5, color=colors[i], label=file_stem)
    
    plt.xlabel('Temperature (T)', fontsize=14)
    plt.ylabel('Heat Capacity (Cv)', fontsize=14)
    plt.title('Combined Heat Capacity vs Temperature', fontsize=16)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'combined_Cv.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create combined C1 plot
    plt.figure(figsize=(10, 7))
    for i, (file_stem, data) in enumerate(all_data.items()):
        T = data["Temperature"].values
        C1 = data["C1/N"].values
        plt.plot(T, C1, linestyle='-', marker=markers[i % len(markers)], 
                markersize=5, linewidth=1.5, color=colors[i], label=file_stem)
    
    plt.xlabel('Temperature (T)', fontsize=14)
    plt.ylabel('Fraction of connected monomers (chain 1)', fontsize=14)
    plt.title('Combined Chain 1 Contacts vs Temperature', fontsize=16)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'combined_C1.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create combined C2 plot
    plt.figure(figsize=(10, 7))
    for i, (file_stem, data) in enumerate(all_data.items()):
        T = data["Temperature"].values
        C2 = data["C2/N"].values
        plt.plot(T, C2, linestyle='-', marker=markers[i % len(markers)], 
                markersize=5, linewidth=1.5, color=colors[i], label=file_stem)
    
    plt.xlabel('Temperature (T)', fontsize=14)
    plt.ylabel('Fraction of connected monomers (chain 2)', fontsize=14)
    plt.title('Combined Chain 2 Contacts vs Temperature', fontsize=16)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'combined_C2.png'), dpi=300, bbox_inches='tight')
    plt.close()
    
    # Create combined C3 plot
    plt.figure(figsize=(10, 7))
    for i, (file_stem, data) in enumerate(all_data.items()):
        T = data["Temperature"].values
        C3 = data["C3/N"].values
        plt.plot(T, C3, linestyle='-', marker=markers[i % len(markers)], 
                markersize=4, linewidth=1, color=colors[i], label=file_stem)
    
    plt.xlabel('Temperature (T)', fontsize=14)
    plt.ylabel('Fraction of inter-chain connections', fontsize=14)
    plt.title('Combined Inter-Chain Contacts vs Temperature', fontsize=16)
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'combined_C3.png'), dpi=300, bbox_inches='tight')
    plt.close()

def process_file(file_path, output_dir):
    # Extract filename without extension for naming plots
    file_stem = Path(file_path).stem
    data = pd.read_csv(file_path)
    
    # Extract all relevant columns
    T = data["Temperature"].values
    E = data["E"].values
    Esquared = data["Esquared"].values
    Binder = data["Binder"].values
    C1 = data["C1"].values
    C1fluc = data["C1fluc"].values
    C2 = data["C2"].values
    C2fluc = data["C2fluc"].values
    C3 = data["C3"].values
    C3fluc = data["C3fluc"].values
    Cv = data["Cv"].values
    
    
    # Create subdirectory for this file's plots
    file_output_dir = os.path.join(output_dir, file_stem)
    os.makedirs(file_output_dir, exist_ok=True)
    
    # Generate and save plots for all variables
    plot_data(T, E, "System Energy", os.path.join(file_output_dir, "E_vs_T.png"), 
              T_min=0.01, T_max=2.0, T_step=0.01, is_d_plot=True)
    plot_data(T, Esquared, "Energy Squared", os.path.join(file_output_dir, "Esquared_vs_T.png"), 
              T_min=0.01, T_max=2.0, T_step=0.01, is_d_plot=False)
    plot_data(T, Binder, "Binder Cumulant", os.path.join(file_output_dir, "Binder_vs_T.png"), 
              T_min=0.01, T_max=2.0, T_step=0.01, is_d_plot=False)
    plot_data(T, C1, "Fraction of connected monomers (chain 1)", os.path.join(file_output_dir, "C1_vs_T.png"), 
              T_min=0.01, T_max=2.0, T_step=0.01, is_d_plot=True)
    plot_data(T, C1fluc, "Chain 1-Wall Contact Fluctuations", os.path.join(file_output_dir, "C1fluc_vs_T.png"), 
              T_min=0.01, T_max=2.0, T_step=0.01, is_d_plot=False)
    plot_data(T, C2, "Fraction of connected monomers (chain 2)", os.path.join(file_output_dir, "C2_vs_T.png"), 
              T_min=0.01, T_max=2.0, T_step=0.01, is_d_plot=True)
    plot_data(T, C2fluc, "Chain 2-Wall Contact Fluctuations", os.path.join(file_output_dir, "C2fluc_vs_T.png"), 
              T_min=0.01, T_max=2.0, T_step=0.01, is_d_plot=False)
    plot_data(T, C3, "Fraction of inter-chain connections", os.path.join(file_output_dir, "C3_vs_T.png"), 
              T_min=0.01, T_max=2.0, T_step=0.01, is_d_plot=True)
    plot_data(T, C3fluc, "Inter-Chain Contact Fluctuations", os.path.join(file_output_dir, "C3fluc_vs_T.png"), 
              T_min=0.01, T_max=2.0, T_step=0.01, is_d_plot=False)
    plot_data(T, Cv, "Heat Capacity", os.path.join(file_output_dir, "Cv_vs_T.png"), 
              T_min=0.01, T_max=2.0, T_step=0.01, is_d_plot=True)
    
    
    return file_stem, data

def main():
    # Directory containing the result files
    input_dir = '/home/amitkumar/Documents/mpconcmz/build/quadraticpotential'
    # Directory to save all plots
    output_dir = '/home/amitkumar/Documents/mpconcmz/build/plotscobines'
    
    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all CSV files in the input directory
    csv_files = [f for f in os.listdir(input_dir) if f.endswith('.csv')]
    
    if not csv_files:
        print(f"No CSV files found in {input_dir}")
        return
    
    print(f"Found {len(csv_files)} CSV files to process")
    
    # Dictionary to store all data for combined plots
    all_data = {}
    
    # Process each file
    for i, csv_file in enumerate(csv_files, 1):
        file_path = os.path.join(input_dir, csv_file)
        print(f"Processing file {i}/{len(csv_files)}: {csv_file}")
        try:
            file_stem, data = process_file(file_path, output_dir)
            all_data[file_stem] = data
            print(f"Successfully processed {csv_file}")
        except Exception as e:
            print(f"Error processing {csv_file}: {str(e)}")
    
    # Create combined plots if we have multiple files
    if len(all_data) > 1:
        print("\nCreating combined plots...")
        create_combined_plots(all_data, output_dir)
    
    print("\nAll files processed!")

if __name__ == "__main__":
    main()