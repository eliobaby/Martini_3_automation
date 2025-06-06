import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Check that the CSV exists
csv_path = "delta_testing_653.csv"
# Read the CSV
df = pd.read_csv(csv_path)

# Extract columns
x_minh = df["delta_G_Minh"].values
x_tristan = df["delta_G_tristan"].values
y_exp = df["delta_G_exp"].values

def plot_with_regression(x, y, xlabel, ylabel, title, output_filename):
    # Determine common axis limits
    vmin = min(x.min(), y.min())
    vmax = max(x.max(), y.max())

    # Compute linear regression (slope and intercept)
    slope, intercept = np.polyfit(x, y, 1)
    y_pred = slope * x + intercept

    # Compute R^2
    ss_res = np.sum((y - y_pred) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    r2 = 1 - (ss_res / ss_tot) if ss_tot != 0 else np.nan

    # Plot
    fig, ax = plt.subplots(figsize=(6, 6), dpi=1000)
    ax.scatter(x, y, s=5, label="Data points")
    # Regression line
    line_x = np.array([vmin, vmax])
    line_y = slope * line_x + intercept
    ax.plot(line_x, line_y, color='red', linewidth=1, label=f"Fit: y={slope:.2f}x+{intercept:.2f}")

    ax.set_xlim(vmin, vmax)
    ax.set_ylim(vmin, vmax)
    ax.set_aspect('equal', 'box')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{title}\n$R^2$ = {r2:.3f}")
    ax.legend(loc="upper left", fontsize="small")

    fig.savefig(output_filename, dpi=1000)
    plt.close(fig)

# Plot and save for delta_G_Minh vs delta_G_exp
plot_with_regression(
    x_minh,
    y_exp,
    xlabel="delta_G_Minh",
    ylabel="delta_G_exp",
    title="δG MinH vs δG Exp",
    output_filename="MINH.png"
)

# Plot and save for delta_G_tristan vs delta_G_exp
plot_with_regression(
    x_tristan,
    y_exp,
    xlabel="delta_G_tristan",
    ylabel="delta_G_exp",
    title="δG Tristan vs δG Exp",
    output_filename="TRISTAN.png"
)

print("Saved MINH.png and TRISTAN.png with regression lines and R² annotations.")
