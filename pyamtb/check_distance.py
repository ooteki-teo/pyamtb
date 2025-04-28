from .tight_binding_model import read_poscar
import numpy as np
import os
import tomlkit

def calculate_distances(poscar_filename, element1="Mn", element2="O"):
    """
    计算两种元素之间的所有距离
    
    参数:
        poscar_filename (str): POSCAR文件路径
        element1 (str): 第一种元素类型
        element2 (str): 第二种元素类型
        
    返回:
        list: 包含距离信息的列表，按距离排序
    """
    # 读取POSCAR文件
    poscar_data = read_poscar(poscar_filename)
    
    # 获取晶格向量
    lattice = poscar_data['lattice']
    
    
    # 获取指定元素类型的原子索引和坐标
    atom_indices1 = [i for i, element in enumerate(poscar_data["atom_symbols"]) if element == element1]
    atom_indices2 = [i for i, element in enumerate(poscar_data["atom_symbols"]) if element == element2]
    
    # 存储所有距离信息
    all_distances = []
    
    # 计算所有可能的原子对之间的距离
    for idx1, atom1_index in enumerate(atom_indices1):
        coord1 = poscar_data["coordinates"][atom1_index]
        
        for idx2, atom2_index in enumerate(atom_indices2):
            coord2 = poscar_data["coordinates"][atom2_index]
            
            # 考虑周期性边界条件
            # 生成可能的格矢量偏移
            for dx in range(-1, 2):
                for dy in range(-1, 2):
                    for dz in range(-1, 2):
                        R = np.array([dx, dy, dz])
                        
                        # 计算考虑周期性的距离
                        diff = coord2 + R - coord1
                        
                        # 转换为笛卡尔坐标计算实际距离
                        diff_cart = np.dot(diff, lattice)
                        distance = np.linalg.norm(diff_cart)
                        
                        all_distances.append({
                            "atom1_index": atom1_index,
                            "atom2_index": atom2_index,
                            "element1": element1,
                            "element2": element2,
                            "distance": distance,
                            "R_vector": R
                        })
    
    # 按距离排序
    all_distances.sort(key=lambda x: x["distance"])
    
    # 按距离统计
    distance_stats = {}
    for dist_info in all_distances:
        # 四舍五入到小数点后4位，便于统计
        rounded_dist = round(dist_info["distance"], 4)
        if rounded_dist not in distance_stats:
            distance_stats[rounded_dist] = {
                "count": 0,
                "pairs": []
            }
        distance_stats[rounded_dist]["count"] += 1
        distance_stats[rounded_dist]["pairs"].append((dist_info["atom1_index"], dist_info["atom2_index"], dist_info["R_vector"]))
    
    # 将统计结果按距离排序
    sorted_stats = sorted(distance_stats.items())
    print("\n距离统计结果:")
    print(f"{'距离(a.u.)':<12}{'出现次数':<10}{'原子对示例'}")
    print("-" * 50)
    for dist, info in sorted_stats:
        example_pair = info["pairs"][0]
        print(f"{dist:<12.4f}{info['count']:<10}{example_pair}")
    
    return all_distances

def print_distances(distances, max_count=100):
    """
    打印距离信息
    
    参数:
        distances (list): 距离信息列表
        max_count (int): 最多显示的距离数量
    """
    print(f"{'原子1':<8}{'原子2':<8}{'元素1':<6}{'元素2':<6}{'距离(Å)':<10}{'距离(a.u.)':<10}{'格矢量':<15}")
    print("-" * 70)
    
    for i, dist_info in enumerate(distances[:max_count]):
        print(f"{dist_info['atom1_index']:<8}{dist_info['atom2_index']:<8}"
              f"{dist_info['element1']:<6}{dist_info['element2']:<6}"
              f"{dist_info['distance']:<10.4f}{dist_info['distance']:<10.4f}"
              f"{dist_info['R_vector']}")

if __name__ == "__main__":
    # 示例用法
    # structure_dir = "structure_to_model"
    # for filename in os.listdir(structure_dir):
    #     if filename.endswith(".vasp"):
    #         print(f"\n处理文件: {filename}")
    poscar_path = os.path.join(".", "kagome_breathing.vasp")
    distances = calculate_distances(poscar_path, "Cu", "Cu")
    # print_distances(distances, max_count=40)
