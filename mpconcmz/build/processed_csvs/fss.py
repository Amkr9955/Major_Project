import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# ======= CONFIG ========
file_configs = [
    {
        "label": "N=10",
        "file_path": "/home/amitkumar/Documents/mpconcmz/build/processed_csvs/results10quadraticpotentialcc1cw1.csv",
        "droplines": [(0.27, "Tm"), (0.675, "Tdes")],
        "marker": "o",   # hollow circle
        "color": "black"
    },
    {
        "label": "N=30",
        "file_path": "/home/amitkumar/Documents/mpconcmz/build/processed_csvs/results30quadraticpotentialcc1cw1.csv",
        "droplines": [(0.34, "Tm"), (0.77, "Tdes")],
        "marker": "s",   # hollow square
        "color": "Steelblue"
    }
]

# Define droplines for each observable and system
droplines_config = {
    "Cv": {
        "N=10": [(0.27, "Tm"), (0.675, "Tdes")],
        "N=30": [(0.34, "Tm"), (0.77, "Tdes")]
    },
    "C1": {
        "N=10": [(0.675, "Tdes")],
        "N=30": [(0.77, "Tdes")]
    },
    "C2": {
        "N=10": [(0.675, "Tdes")],
        "N=30": [(0.77, "Tdes")]
    },
    "C3": {
        "N=10": [(0.27, "Tm")],
        "N=30": [(0.34, "Tm")]
    },
    "E": {
        "N=10": [(0.27, "Tm"), (0.675, "Tdes")],
        "N=30": [(0.34, "Tm"), (0.77, "Tdes")]
    }
}

T_min, T_max, T_step = 0.01, 2.0, 0.01

# List of observables to plot
observables = [
    ("E", "<E>/N", True),
    ("Cv", "Cv", True),
    ("C1", "<Nc>/N (chain 1)", True),
    ("C2", "<Nc>/N (chain 2)", True),
    ("C3", "<Nc>/N (chain 1-chain 2)", True),
]

def plot_comparison(observable_key, y_label, filename, is_d_plot=True):
    plt.figure(figsize=(10, 7) if is_d_plot else (7, 5))

    y_min_total, y_max_total = float('inf'), -float('inf')

    for cfg in file_configs:
        data = pd.read_csv(cfg["file_path"])
        T = data["Temperature"].values
        Y = data[observable_key].values

        mask = (T >= T_min) & (T <= T_max)
        T_filtered, Y_filtered = T[mask], Y[mask]

        if T_step is not None and len(T_filtered) > 1:
            unique_T = np.unique(T_filtered)
            step = int(T_step / np.min(np.diff(unique_T)))
            selected_T = unique_T[::step] if step > 0 else unique_T
            step_mask = np.isin(T_filtered, selected_T)
            T_filtered, Y_filtered = T_filtered[step_mask], Y_filtered[step_mask]

        sorted_indices = np.argsort(T_filtered)
        T_sorted = T_filtered[sorted_indices]
        Y_sorted = Y_filtered[sorted_indices]

        y_min_total = min(y_min_total, np.min(Y_sorted))
        y_max_total = max(y_max_total, np.max(Y_sorted))

        # Plot curve with hollow markers
        plt.plot(T_sorted, Y_sorted,
                 linestyle='-', marker=cfg["marker"],
                 markerfacecolor='none', markeredgewidth=1.5,
                 markersize=6, color=cfg["color"], linewidth=2,
                 label=cfg["label"])

        # Plot droplines and annotate
                # Get droplines for this observable and system
        system_label = cfg["label"]
        dropline_Ts = droplines_config.get(observable_key, {}).get(system_label, [])

        for T_drop, drop_label in dropline_Ts:

            if T_min <= T_drop <= T_max:
                Y_interp = np.interp(T_drop, T_sorted, Y_sorted)
                plt.plot([T_drop, T_drop], [y_min_total - 0.01, Y_interp],
                         color='green', linestyle='--', linewidth=2.5)
                plt.annotate(drop_label, xy=(T_drop, y_min_total - 0.01),
                             xytext=(0, -12), textcoords='offset points',
                             ha='center', va='top', fontsize=15, color='green')

        # Add annotation box
                     # Add annotation box per system
        if dropline_Ts:
            box_lines = [f"N = {system_label.split('=')[1].strip()}"]  # Extract N value from label like 'N=100'
            for T_val, label in dropline_Ts:
                box_lines.append(f"{label:<4} = {T_val:.3f}")
            box_text = '\n'.join(box_lines)

            # Automatically shift boxes vertically to avoid overlap
            if system_label == "N=10":
                x_pos, y_pos = 0.98, 0.88
            elif system_label == "N=30":
                x_pos, y_pos = 0.98, 0.74
            else:
                x_pos, y_pos = 0.98, 0.74  # Default position

            plt.gca().text(
                x_pos, y_pos, box_text,
                transform=plt.gca().transAxes,
                ha='right', va='top',
                fontsize=11, color='black',
                bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.4')
            )



    # Final styling
    y_padding = 0.05 * (y_max_total - y_min_total)
    plt.ylim(y_min_total - y_padding, y_max_total + y_padding)
    plt.xlabel('Temperature (T)', fontsize=14)
    plt.ylabel(y_label, fontsize=14)
    plt.title(f'{y_label} vs Temperature', fontsize=16)
    plt.grid(True, linestyle='--', linewidth=0.5)
    plt.legend(fontsize=11)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.tight_layout()
    plt.savefig(filename, dpi=300)
    plt.close()

# ======= Run Plots for all observables ========
for key, label, is_d_plot in observables:
    plot_comparison(key, label, f"{key}_compare.png", is_d_plot)

print("All comparative plots saved!")
