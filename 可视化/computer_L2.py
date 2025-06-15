import numpy as np
from io_utils import parse_out_file, compute_L2_error
from plot_utils import plot_L2_convergence

# 文件列表与对应单元尺寸
filenames = [
    r"data\conv1.out",
    r"data\conv4.out",
    r"data\conv16.out"
]
h_values = [np.sqrt(5)/16, np.sqrt(5)/4, np.sqrt(5)]

L2_errors = []

for filename in filenames:
    print(f"Processing {filename}")
    element_type, node_coords, displacements = parse_out_file(filename)
    error = compute_L2_error(element_type, node_coords, displacements)
    print(f"L2 error: {error:.6e}")
    L2_errors.append(error)

# 绘图与计算收敛率
slope = plot_L2_convergence(h_values, L2_errors)
