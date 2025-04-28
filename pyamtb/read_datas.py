import tomlkit
import os
from typing import Optional
import numpy as np

def read_poscar(filename="POSCAR", selected_elements: Optional[list]=None, to_direct: bool=True, use_scale: bool=True):
    """
    读取VASP格式的POSCAR文件，并将结构信息存储在字典中返回
    
    参数:
        filename (str): POSCAR文件路径，默认为"POSCAR"
        selected_elements (list): 可选参数，用于指定要读取的元素类型列表
        to_direct (bool): 可选参数，用于指定是否将坐标转换为分数坐标
        use_scale (bool): 可选参数，用于指定是否使用比例因子
    返回:
        dict: 包含结构信息的字典，包括以下键:
            - comment: 注释行
            - scale: 比例因子
            - lattice: 晶格向量 (3x3 numpy数组)
            - elements: 元素类型列表
            - atom_counts: 各元素原子数量列表
            - total_atoms: 总原子数
            - is_direct: 是否为分数坐标 (布尔值)
            - coords: 原子坐标数组 (nx3 numpy数组)
            - atom_symbols: 每个原子的元素符号列表
    """
    try:
        with open(filename, 'r', encoding='utf-8') as f:
            lines = f.readlines()
            if not lines:
                raise ValueError(f"文件 {filename} 为空")

        # 解析基本信息
        comment = lines[0].strip()
        scale = float(lines[1].strip()) if use_scale else 1
        lattice = np.array([list(map(float, line.split())) for line in lines[2:5]]) * scale
        scale = 1
        elements = lines[5].split()
        atom_counts = list(map(int, lines[6].split()))
        total_atoms = sum(atom_counts)
        
        # 判断坐标类型
        coord_type = lines[7].strip()
        if coord_type.lower().startswith('d'):
            is_direct = True
        elif coord_type.lower().startswith('c') or coord_type.lower().startswith('k'):
            is_direct = False
        else:
            raise ValueError(f"未知的坐标类型: {coord_type}")

        # 读取原子坐标
        coords = np.array([list(map(float, line.split()[:3])) for line in lines[8:8+total_atoms]])

        # 如果是笛卡尔坐标，转换为分数坐标
        if not is_direct and to_direct:
            change_to_direct = True
            inv_lattice = np.linalg.inv(lattice)
            coords = np.dot(coords, inv_lattice)
        # 如果是分数坐标，转换为笛卡尔坐标
        elif is_direct and not to_direct:
            change_to_direct = False
            coords = np.dot(coords, lattice)
            is_direct = False
        else:
            change_to_direct = "Do not change"

        # 生成每个原子的元素符号列表
        atom_symbols = []
        for symbol, count in zip(elements, atom_counts):
            atom_symbols.extend([symbol] * count)
        
        # 如果指定了选定元素，则只保留这些元素对应的原子
        if selected_elements is not None:
            selected_indices = []
            selected_atom_symbols = []
            selected_coords = []
            
            for i, symbol in enumerate(atom_symbols):
                if symbol in selected_elements:
                    selected_indices.append(i)
                    selected_atom_symbols.append(symbol)
                    selected_coords.append(coords[i])
            
            # 更新原子信息
            atom_symbols = selected_atom_symbols
            coords = np.array(selected_coords) if selected_coords else np.empty((0, 3))
            total_atoms = len(selected_indices)
            
            # 更新元素计数
            new_elements = []
            new_atom_counts = []
            for element in selected_elements:
                if element in elements:
                    count = atom_symbols.count(element)
                    if count > 0:
                        new_elements.append(element)
                        new_atom_counts.append(count)
            
            elements = new_elements
            atom_counts = new_atom_counts

        # 构建结果字典
        poscar_data = {
            "comment": comment,
            "scale": scale,
            "lattice": lattice,
            "elements": elements,
            "atom_counts": atom_counts,
            "total_atoms": total_atoms,
            "is_direct": is_direct,
            "coordinates": coords,
            "atom_symbols": atom_symbols
        }
        # print(poscar_data)
        return poscar_data
    except Exception as e:
        print(f"读取POSCAR文件时出错: {str(e)}")
        print(f"文件内容:")
        with open(filename, 'r', encoding='utf-8') as f:
            print(f.read())
        raise


# 计算字符串的字符频率向量
def get_char_freq(s):
    freq = [0] * 26
    s = s.lower()
    for c in s:
        if c.isalpha():
            freq[ord(c) - ord('a')] += 1
    return np.array(freq)

# 计算余弦相似度
def cosine_similarity(v1, v2):
    dot = np.dot(v1, v2)
    norm1 = np.linalg.norm(v1)
    norm2 = np.linalg.norm(v2)
    if norm1 == 0 or norm2 == 0:
        return 0
    return dot / (norm1 * norm2)

def find_similar_strings(target: str, string_list: list) -> list:
    """
    从字符串列表中找到与目标字符串相似的字符串
    
    参数:
        target (str): 目标字符串
        string_list (list): 字符串列表
        
    返回:
        list: 相似字符串列表
    """
    # 将目标字符串转为小写以进行不区分大小写的比较
    target = target.lower()
    
    # 计算目标字符串的频率向量
    target_freq = get_char_freq(target)
    
    # 找到相似的字符串
    similar_strings = []
    for s in string_list:
        s_freq = get_char_freq(s)
        similarity = cosine_similarity(target_freq, s_freq)
        print(f"相似度: {similarity}, 字符串: {s} -> {target}")
        if similarity > 0.8:  # 设置相似度阈值
            similar_strings.append(s)
    
    return similar_strings


def read_parameters(filename="tbparas.toml"):
    """
    从toml文件中读取参数
    
    参数:
        filename (str): toml文件路径
        
    返回:
        dict: 包含参数的字典
    """
    if not os.path.exists(filename):
        raise FileNotFoundError(f"找不到参数文件: {filename}")
        
    with open(filename, 'r', encoding='utf-8') as f:
        params = tomlkit.load(f)
    
    # 设置默认值
    default_params = {
        "poscar_filename": "Mn2N.vasp",
        "lattice_constant": 1.0,
        "t0": 1.0,
        "t0_distance": 2.0,
        "min_distance": 0.1,
        "max_distance": 2.6,
        "dimk": 2,
        "dimr": 3,
        "hopping_decay": 1,
        "use_elements": ["Mn", "N"],
        "output_filename": "Mn2N_model.png",
        "output_format": "png",
        "savedir": ".",
        "same_atom_negative_coupling": False,
        "magnetic_moment": 0.1,
        "magnetic_order": "+-0",
        "nspin": 2,
        "onsite_energy": [0.0, 0.0, 0.0],
        "ylim": [-1, 1],
        "kpath": [[0, 0], [0.5, 0], [0.5, 0.5], [0.0, 0.5], [0.0, 0.0], [0.5, 0.5]],
        "klabel": ["G", "X", "M", "Y", "G", "M"],
        "nkpt": 100,
        "max_neighbors": 1,
        "is_print_tb_model_hop": True,
        "is_check_flat_bands": True,
        "is_print_tb_model": True,
        "is_black_degenerate_bands": True,
        "energy_threshold": 1e-5
    }
    
    # 用文件中的值更新默认值
    for key in default_params:
        if key in params:
            default_params[key] = params[key]
    
    # 检查是否有未知参数
    for key in params.keys():
        # print(f"key: {key}")
        if key not in default_params.keys():
            # 查找相似的参数名
            similar_keys = find_similar_strings(key, default_params.keys())
            if similar_keys:
                raise ValueError(f"输入参数{key}有误，是否想使用以下参数之一: {similar_keys}?")
            else:
                raise ValueError(f"输入参数{key}有误，未找到类似参数，请仔细检查")
    
    # 几个需要为整数的参数
    default_params["dimk"] = int(default_params["dimk"])
    default_params["dimr"] = int(default_params["dimr"])
    default_params["nkpt"] = int(default_params["nkpt"])
    default_params["max_neighbors"] = int(default_params["max_neighbors"])

    poscar=read_poscar(default_params["poscar_filename"])
    for ele in default_params["use_elements"]:
        if ele not in poscar["elements"]:
            raise ValueError(f"{filename} 的 use_elements 参数有误，元素{ele}不在POSCAR中")

    return default_params
