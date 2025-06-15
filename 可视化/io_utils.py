import re
import numpy as np

def parse_out_file(filename, debug=0):
    with open(filename, 'r') as f:
        lines = f.readlines()

    node_start = elem_start = disp_start = None
    for i, line in enumerate(lines):
        if 'N O D A L   P O I N T   D A T A' in line:
            node_start = i+3
        if 'E L E M E N T   I N F O R M A T I O N' in line:
            elem_start = i+2
        if 'D I S P L A C E M E N T S' in line:
            disp_start = i+2
        if 'ELEMENT TYPE' in line:
            m = re.search(r'=\s*(\d+)', line)
            element_type = int(m.group(1))

    node_coords = read_nodes(lines, node_start, element_type)
    displacements = read_displacements(lines, disp_start, element_type)

    return element_type, node_coords, displacements

def read_nodes(lines, start, element_type):
    node_coords = []
    for line in lines[start:]:
        if re.match(r'\s*\d+\s+\d+\s+\d+', line):
            parts = line.split()
            coords = list(map(float, parts[-3:] if element_type == 1 else parts[-2:]))
            node_coords.append(coords)
        elif line.strip() == '':
            break
    return np.array(node_coords)

def read_displacements(lines, start, element_type):
    displacements = []
    for line in lines[start:]:
        if re.match(r'\s*\d+\s+[-\deE.+]+\s+[-\deE.+]+', line):
            parts = line.split()
            disp = list(map(float, parts[1:4])) if element_type == 1 else list(map(float, parts[1:3]))
            displacements.append(disp)
        elif line.strip() == '':
            break
    return np.array(displacements)

def compute_L2_error(element_type, node_coords, displacements):
    if element_type == 2:
        u_exact = 0.12 * node_coords[:, 0] * (node_coords[:, 1] - 0.5)
        v_exact = -0.06 * node_coords[:, 0]**2
        u_num = displacements[:, 0]
        v_num = displacements[:, 1]
        error = np.sqrt(np.mean((u_num - u_exact)**2 + (v_num - v_exact)**2))
        return error
    else:
        raise ValueError("Only Q4 (element_type=2) supported for L2 computation.")
