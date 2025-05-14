import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Load the data
file_path = '/home/amitkumar/Documents/mpconcmz/build/processed_csvs/results20quadraticpotentialcc1cw1sideways.csv'
data = pd.read_csv(file_path)

# Extract all relevant columns
T = data["Temperature"].values
E = data["E"].values
Esquared = data["Esquared"].values
Binder =data["Binder"].values
C1 = data["C1"].values
C1fluc = data["C1fluc"].values
C2 = data["C2"].values
C2fluc = data["C2fluc"].values
C3 = data["C3"].values
C3fluc = data["C3fluc"].values
Cv = data["Cv"].values

def plot_data(T, Y, y_label, filename, T_min=None, T_max=None, T_step=None, is_d_plot=False, dropline_Ts=None):
    import matplotlib.ticker as ticker

    # Create mask for temperature range
    mask = np.ones_like(T, dtype=bool)
    if T_min is not None:
        mask = mask & (T >= T_min)
    if T_max is not None:
        mask = mask & (T <= T_max)

    T_filtered = T[mask]
    Y_filtered = Y[mask]

    if T_step is not None and len(T_filtered) > 0:
        unique_T = np.unique(T_filtered)
        step = int(T_step / np.min(np.diff(unique_T)))
        selected_T = unique_T[::step] if step > 0 else unique_T
        step_mask = np.isin(T_filtered, selected_T)
        T_filtered = T_filtered[step_mask]
        Y_filtered = Y_filtered[step_mask]

    sorted_indices = np.argsort(T_filtered)
    T_sorted, Y_sorted = T_filtered[sorted_indices], Y_filtered[sorted_indices]

    # Determine dynamic y-axis limits
    y_min_data = np.min(Y_sorted)
    y_max_data = np.max(Y_sorted)
    y_padding = 0.05 * (y_max_data - y_min_data)
    y_bottom = y_min_data - y_padding  # Just a little below the lowest point
    y_top = y_max_data + y_padding     # Just a little above the highest point

    plt.figure(figsize=(10, 7) if is_d_plot else (7, 5))
    plt.plot(T_sorted, Y_sorted, color='black', linestyle="-", marker='o',
             markersize=4 if is_d_plot else 3, linewidth=2 if is_d_plot else 1, label='Data')

    plt.xlabel('Temperature (T)', fontsize=14 if is_d_plot else 12)
    plt.ylabel(y_label, fontsize=14 if is_d_plot else 12)
    plt.title(f'{y_label} vs Temperature', fontsize=16 if is_d_plot else 12)
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.xticks(fontsize=12 if is_d_plot else 10)
    plt.yticks(fontsize=12 if is_d_plot else 10)
    plt.legend(fontsize=12 if is_d_plot else 10)
    plt.ylim(y_bottom, y_top)

    # Dropline logic: tuple of (T, label)
    if dropline_Ts:
        # Plot droplines with separate x-axis labels
        for T_drop, label in dropline_Ts:
            if T_drop < T_sorted[0] or T_drop > T_sorted[-1]:
                continue
            Y_interp = np.interp(T_drop, T_sorted, Y_sorted)
            plt.plot([T_drop, T_drop], [y_bottom, Y_interp], color='green', linestyle='--', linewidth=2.5)

            # Label below each line
            plt.annotate(label, xy=(T_drop, y_bottom), xytext=(0, -14), textcoords='offset points',
                         ha='center', va='top', fontsize=15, color='green', rotation=0)

        # ✅ Construct annotation box text
        box_lines = ["N=20"]  # Add observable name
        for T_val, label in dropline_Ts:
            box_lines.append(f"{label:<4} = {T_val:.3f}")  # Fixed width for alignment

        box_text = '\n'.join(box_lines)

        # ✅ Place the box (top-right, black text)
        plt.gca().text(
            0.98, 0.92, box_text,
            transform=plt.gca().transAxes,
            ha='right', va='top',
            fontsize=11, color='black',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.4')
        )


    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()



droplines1 = [(0.06,"Tads"),(0.34,"Tm"), (0.70,"Tdes")] 
droplines2 = [(0.34,"Tm")]
droplines3 = [(0.06,"Tads"),(0.70,"Tdes")]
droplines4 = [(0.70,"Tdes")]


# Generate and save plots for all variables
plot_data(T, E, "<E>", "E_vs_T.png", T_min=0.01, T_max=2.0, T_step=0.01, is_d_plot=True,dropline_Ts=droplines1)
plot_data(T, C1, "<Nc>/N(chain 1)", "C1_vs_T.png", T_min=0.01, T_max=2.0, T_step=0.01, is_d_plot=True,dropline_Ts=droplines4)
plot_data(T, C2, "<Nc>/N(chain 2)", "C2_vs_T.png", T_min=0.01, T_max=2.0, T_step=0.01, is_d_plot=True,dropline_Ts=droplines3)
plot_data(T, C3, "<Nc>/N (inter-chain)", "C3_vs_T.png", T_min=0.01, T_max=2.0, T_step=0.01, is_d_plot=True,dropline_Ts=droplines2)
plot_data(T, Cv, "Cv", "Cv_vs_T.png", T_min=0.01, T_max=2.0, T_step=0.01, is_d_plot=True,dropline_Ts=droplines1)


print("All plots saved successfully!")