import re
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from pathlib import Path

# --- 配置区 (Configuration Area) ---
# 将调试开关和默认数据目录放在这里，方便修改
DEBUG = False
DEFAULT_DATA_DIR = Path("./data")  # 使用相对路径，或改为 Path("E:/Element/src/data")


# --- 定义抽象基类 (Define Abstract Base Class) ---
class FEADataExtractor:
    """
    有限元分析结果提取器的抽象基类。
    定义了所有单元类型处理器都应遵循的通用接口和共享逻辑。
    """
    # 定义了 .out 文件中各个数据块的起始标记
    SECTION_MARKERS = {
        'nodes': 'N O D A L   P O I N T   D A T A',
        'elements': 'E L E M E N T   I N F O R M A T I O N',
        'displacements': 'D I S P L A C E M E N T S',
        'element_type': 'ELEMENT TYPE'
    }

    def __init__(self, lines):
        """
        构造函数。
        :param lines: 从 .out 文件读取的所有行。
        """
        self.lines = lines
        self.section_starts = self._find_section_starts()
        # 初始化属性
        self.is_3d = False
        self.nodes_per_element = 0

    def _find_section_starts(self):
        """在文件中查找所有关键数据块的起始行号。"""
        starts = {}
        for i, line in enumerate(self.lines):
            for key, marker in self.SECTION_MARKERS.items():
                if marker in line:
                    # 根据不同部分的头文件行数进行偏移
                    offset = 3 if key == 'nodes' else 2
                    starts[key] = i + offset
        return starts

    def _parse_data_block(self, start_line, parser_func):
        """通用的数据块解析逻辑。"""
        data = []
        if start_line is None:
            return data
        for line in self.lines[start_line:]:
            if not line.strip():  # 遇到空行则停止
                break
            match = parser_func(line)
            if match:
                data.append(match)
        return data

    def get_element_type(self):
        """从文件中解析单元类型编号。"""
        line_index = self.section_starts.get('element_type', 0) - 2
        if line_index > 0:
            match = re.search(r'=\s*(\d+)', self.lines[line_index])
            if match:
                return int(match.group(1))
        return None

    # --- 以下为需要子类实现的抽象方法 ---
    def read_nodes(self):
        raise NotImplementedError("子类必须实现 read_nodes 方法。")

    def read_elements(self):
        raise NotImplementedError("子类必须实现 read_elements 方法。")

    def read_displacements(self):
        raise NotImplementedError("子类必须实现 read_displacements 方法。")

    def plot_mesh(self, ax, coords, elements, title, color):
        raise NotImplementedError("子类必须实现 plot_mesh 方法。")


# --- 定义具体单元类型的子类 (Define Concrete Subclasses) ---
class TrussExtractor(FEADataExtractor):
    """桁架单元 (Truss) 处理器。"""

    def __init__(self, lines):
        super().__init__(lines)
        self.is_3d = True
        self.nodes_per_element = 2

    def read_nodes(self):
        """读取节点坐标 (x, y, z)。"""

        def parser(line):
            parts = line.split()
            if len(parts) >= 4 and re.match(r'\s*\d+', line):
                return [float(p) for p in parts[-3:]]
            return None

        coords = self._parse_data_block(self.section_starts.get('nodes'), parser)
        return np.array(coords)

    def read_elements(self):
        """读取单元连接关系 (node1, node2)。"""

        def parser(line):
            parts = line.split()
            if len(parts) >= 4 and re.match(r'\s*\d+', line):
                return [int(p) for p in parts[1:self.nodes_per_element + 1]]
            return None

        return self._parse_data_block(self.section_starts.get('elements'), parser)

    def read_displacements(self):
        """读取节点位移 (dx, dy, dz)。"""

        def parser(line):
            parts = line.split()
            if len(parts) >= 4 and re.match(r'\s*\d+', line):
                return [float(p) for p in parts[1:4]]
            return None

        displacements = self._parse_data_block(self.section_starts.get('displacements'), parser)
        return np.array(displacements)

    def plot_mesh(self, ax, coords, elements, title, color):
        """绘制三维网格。"""
        ax.set_title(title)
        for elem in elements:
            # 节点编号从1开始，需减1转换为索引
            pts = coords[np.array(elem) - 1]
            ax.plot(pts[:, 0], pts[:, 1], pts[:, 2], f'{color}-')
        ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2], c=color, s=50)  # s调整点大小

        # 设置坐标轴比例
        max_range = np.array([coords[:, 0].max() - coords[:, 0].min(),
                              coords[:, 1].max() - coords[:, 1].min(),
                              coords[:, 2].max() - coords[:, 2].min()]).max() / 2.0
        mid_x = (coords[:, 0].max() + coords[:, 0].min()) * 0.5
        mid_y = (coords[:, 1].max() + coords[:, 1].min()) * 0.5
        mid_z = (coords[:, 2].max() + coords[:, 2].min()) * 0.5
        ax.set_xlim(mid_x - max_range, mid_x + max_range)
        ax.set_ylim(mid_y - max_range, mid_y + max_range)
        ax.set_zlim(mid_z - max_range, mid_z + max_range)


class Q4Extractor(FEADataExtractor):
    """四节点四边形单元 (Q4) 处理器。"""

    def __init__(self, lines):
        super().__init__(lines)
        self.is_3d = False
        self.nodes_per_element = 4

    def read_nodes(self):
        """读取节点坐标 (x, y)。"""

        def parser(line):
            parts = line.split()
            if len(parts) >= 3 and re.match(r'\s*\d+', line):
                return [float(p) for p in parts[-2:]]
            return None

        coords = self._parse_data_block(self.section_starts.get('nodes'), parser)
        return np.array(coords)

    def read_elements(self):
        """读取单元连接关系 (node1, node2, node3, node4)。"""

        def parser(line):
            parts = line.split()
            if len(parts) >= 6 and re.match(r'\s*\d+', line):
                return [int(p) for p in parts[1:self.nodes_per_element + 1]]
            return None

        return self._parse_data_block(self.section_starts.get('elements'), parser)

    def read_displacements(self):
        """读取节点位移 (dx, dy)。"""

        def parser(line):
            parts = line.split()
            if len(parts) >= 3 and re.match(r'\s*\d+', line):
                return [float(p) for p in parts[1:3]]
            return None

        displacements = self._parse_data_block(self.section_starts.get('displacements'), parser)
        return np.array(displacements)

    def plot_mesh(self, ax, coords, elements, title, color):
        """绘制二维网格。"""
        ax.set_title(title)
        ax.axis('equal')
        for elem in elements:
            pts = coords[np.array(elem) - 1]
            pts = np.vstack([pts, pts[0]])  # 闭合图形
            ax.plot(pts[:, 0], pts[:, 1], f'{color}-')
        ax.scatter(coords[:, 0], coords[:, 1], c=color)


# --- 工厂函数 (Factory Function) ---
def get_extractor(element_type, lines):
    """根据单元类型编号返回相应的处理器实例。"""
    if element_type == 1:
        return TrussExtractor(lines)
    elif element_type == 2:
        return Q4Extractor(lines)
    else:
        raise ValueError(f"不支持的单元类型: {element_type}")


# --- 主程序 (Main Execution) ---
def main():
    """主函数，控制程序流程。"""
    try:
        # 1. 获取用户输入
        user_file = input("请输入文件名 (不含路径和扩展名, 如 test): ")
        scale = float(input("请输入位移放大系数 (如 1, 10, 100): "))

        # 使用 pathlib 构建路径，更健壮
        filepath = DEFAULT_DATA_DIR / f"{user_file}.out"
        if not filepath.exists():
            print(f"错误: 文件不存在 -> {filepath}")
            return

        # 2. 读取文件内容
        with open(filepath, 'r') as f:
            lines = f.readlines()

        # 3. 创建合适的处理器
        # 先用一个临时实例来获取单元类型
        temp_extractor = FEADataExtractor(lines)
        element_type = temp_extractor.get_element_type()
        if element_type is None:
            print("错误: 未能在文件中找到 'ELEMENT TYPE'。")
            return

        extractor = get_extractor(element_type, lines)

        # 4. 使用处理器提取数据
        print(f"正在处理单元类型: {extractor.__class__.__name__}")
        node_coords = extractor.read_nodes()
        elements = extractor.read_elements()
        displacements = extractor.read_displacements()

        if DEBUG:
            print("节点坐标:\n", node_coords)
            print("单元连接:\n", elements)
            print("位移:\n", displacements)

        if node_coords.size == 0 or not elements or displacements.size == 0:
            print("错误：数据提取不完整，请检查 .out 文件格式或程序。")
            return

        # 5. 计算变形后坐标
        deformed_coords = node_coords + scale * displacements

        # 6. 绘图
        # 原始结构
        fig1 = plt.figure(figsize=(10, 8))
        ax1 = fig1.add_subplot(111, projection='3d' if extractor.is_3d else None)
        extractor.plot_mesh(ax1, node_coords, elements, 'Original Mesh', color='b')
        img_path1 = filepath.parent / f"{filepath.stem}_original.png"
        plt.savefig(img_path1, dpi=300)
        print(f"原始结构图已保存至: {img_path1}")
        plt.show()

        # 变形后结构
        fig2 = plt.figure(figsize=(10, 8))
        ax2 = fig2.add_subplot(111, projection='3d' if extractor.is_3d else None)
        extractor.plot_mesh(ax2, deformed_coords, elements, f'Deformed Mesh (scale={scale})', color='r')
        img_path2 = filepath.parent / f"{filepath.stem}_deformed.png"
        plt.savefig(img_path2, dpi=300)
        print(f"变形结构图已保存至: {img_path2}")
        plt.show()

    except FileNotFoundError:
        print(f"错误: 找不到输入文件。请确保 '{filepath}' 存在。")
    except ValueError as e:
        print(f"输入错误: {e}")
    except Exception as e:
        print(f"程序运行出错: {e}")


if __name__ == '__main__':
    main()