import pandas as pd
import matplotlib.pyplot as plt
import sys

# Usage: python coverage_plot.py depth.txt output.png

depth_file = sys.argv[1]
output_file = sys.argv[2]

# Load data
df = pd.read_csv(depth_file, sep="\t", header=None)
df.columns = ["chrom", "position", "depth"]

# Plot
plt.figure(figsize=(12, 4))
plt.plot(df["position"], df["depth"])
plt.xlabel("Genome Position")
plt.ylabel("Read Depth")
plt.title("Genome Coverage")
plt.tight_layout()

# Save
plt.savefig(output_file, dpi=300)
plt.close()

print("Coverage plot saved:", output_file)
