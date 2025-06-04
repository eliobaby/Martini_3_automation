import matplotlib.pyplot as plt

# Raw data from the provided GRO snippet
data = """
1TPCN00004
25
    1TPCN0  TC3    1   7.567  -0.801  -3.305
    1TPCN0  TC3    2   6.632  -1.048  -1.345
    1TPCN0 TN4a    3   7.364   0.155   0.281
    1TPCN0  TN6    4   9.435   0.716   0.004
    1TPCN0 TN4a    5  11.584   1.639   0.117
    1TPCN0 TN2a    6  13.393   3.039  -0.673
    1TPCN0  TN6    7  10.984   3.792   0.459
    1TPCN0   N5    8   9.810   2.321   2.124
    1TPCN0  TC5    9   6.298  -1.660   0.693
    1TPCN0  N4a   10   5.288  -4.134   1.034
    1TPCN0  N4a   11   3.076  -1.626  -1.166
    1TPCN0  SC4   12  -0.604  -2.160   0.677
    1TPCN0  SC3   13  -1.228  -0.395  -0.411
    1TPCN0  SP2   14  -4.214  -2.550  -0.144
    1TPCN0  TC3   15  -6.541  -1.263  -0.403
    1TPCN0  TC5   16  -8.205  -1.954   0.780
    1TPCN0  N4a   17  -8.775  -4.447  -0.465
    1TPCN0 TN4a   18  -7.461   0.030   1.045
    1TPCN0  TN6   19  -8.163   1.851   0.149
    1TPCN0 TN4a   20  -9.140   3.297  -0.943
    1TPCN0 TN2a   21  -9.922   2.420  -3.064
    1TPCN0  TN6   22 -11.477   2.550  -0.763
    1TPCN0   N5   23 -10.226   2.729   1.696
    1TPCN0  TC3   24  -5.544   0.710  -1.571
    1TPCN0   N6   25   0.807   1.652   0.884
  10.00000  10.00000  10.00000
"""

# Parse coordinates and labels
lines = data.strip().splitlines()
coords = []
labels = []

for line in lines[2:-1]:
    parts = line.split()
    idx = int(parts[2])
    x, y, z = map(float, parts[3:6])
    coords.append((x, y, z))
    labels.append(str(idx))

# Plotting
fig = plt.figure()
ax = fig.add_subplot(projection='3d')
xs, ys, zs = zip(*coords)
ax.scatter(xs, ys, zs)

# Annotate each point with its index
for x, y, z, label in zip(xs, ys, zs, labels):
    ax.text(x, y, z, label)

ax.set_xlabel('X (Å)')
ax.set_ylabel('Y (Å)')
ax.set_zlabel('Z (Å)')
ax.set_title('3D Structure from BEADS Coordinates')
plt.tight_layout()
plt.show()
