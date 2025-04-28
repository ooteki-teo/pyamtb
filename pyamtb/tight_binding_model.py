"""
读取`tbparas.toml`文件中的参数，利用pythtb库建立紧束缚模型

calculate_distance(coord1, coord2, lattice) 计算两个原子之间的距离,考虑周期性边界条件

hopping_strength(distance) 计算跃迁强度,随距离指数衰减

create_pythtb_model(poscar_filename) 从POSCAR文件创建紧束缚模型,设置轨道位置和跃迁参数

calculate_band_structure(model) 计算能带结构并绘图


"""

from pythtb import * # import TB model class
import numpy as np
import matplotlib.pyplot as plt
import os
from .read_datas import read_poscar
from .parameters import Parameters
from copy import deepcopy

# 创建全局参数实例
params = Parameters()

DEBUG = True

def str_to_mag(mag_str, magnetic_moment):
    """
    将磁序字符串转换为磁矩列表
    """
    return [magnetic_moment if c == '+' else -magnetic_moment if c == '-' else 0 for c in mag_str]

def calculate_distance(coord1, coord2, lattice):
    """
    计算两个原子之间的距离（考虑周期性边界条件）
    
    参数:
        coord1 (np.array): 第一个原子的分数坐标
        coord2 (np.array): 第二个原子的分数坐标
        lattice (np.array): 晶格矩阵
        
    返回:
        float: 两个原子之间的距离（埃）
    """
    # 计算分数坐标差值
    diff = coord2 - coord1
    
    # 应用最小镜像约定（考虑周期性边界条件）
    diff = diff - np.round(diff)
    
    # 转换为笛卡尔坐标
    cart_diff = np.dot(diff, lattice)
    
    # 计算距离
    distance = np.linalg.norm(cart_diff)
    
    return distance

def hopping_strength(distance):
    """
    计算跃迁强度与距离的关系
    
    参数:
        distance (float): 实际距离（埃）
        t0_distance (float): 参考距离（埃）
        t0 (float): 参考距离下的跃迁强度
        
    返回:
        float: 跃迁强度
    """
    # 使用指数衰减模型: t = t0 * exp(-lambda_*(d-d0)/d0)
    return params.t0 * np.exp(-params.lambda_*(distance - params.t0_distance) / params.t0_distance)

def get_coupling_strength(atom1_index, atom2_index, poscar_data):
    """
    计算两个原子之间的耦合强度，考虑周期性边界条件和相邻格点上的等价原子
    
    参数:
        atom1_index (int): 第一个原子的索引
        atom2_index (int): 第二个原子的索引
        poscar_data (dict): POSCAR数据字典
        
    返回:
        tuple: (耦合强度数组, 格矢量数组) - 包含中心格点和相邻格点的耦合信息
    """
    # 获取原子坐标和元素类型
    coord1 = poscar_data["coordinates"][atom1_index]
    coord2 = poscar_data["coordinates"][atom2_index]

    type1 = poscar_data["atom_symbols"][atom1_index]
    type2 = poscar_data["atom_symbols"][atom2_index]

    # 如果两个原子是同一种元素，则耦合强度为正，否则为负（没必要，但是可以）
    coupling_sign = -1 if ((type1 == type2) and params.same_atom_negative_coupling) else 1
    
    # 生成相邻格点的格矢量
    R_vectors = []
    coupling_values = []
    distance_values = []

    # 生成所有可能的格矢量
    Rlist = []
    if params.dimk == 3:
        for i in range(-params.max_neighbors, params.max_neighbors+1):
            for j in range(-params.max_neighbors, params.max_neighbors+1):
                for k in range(-params.max_neighbors, params.max_neighbors+1):
                    Rlist.append(np.array([i, j, k]))
    elif params.dimk == 2:
        for i in range(-params.max_neighbors, params.max_neighbors+1):
            for j in range(-params.max_neighbors, params.max_neighbors+1):
                Rlist.append(np.array([i, j, 0]))
    elif params.dimk == 1:
        Rlist = [[i,0,0] for i in range(-params.max_neighbors, params.max_neighbors+1)]
    else:
        raise ValueError(f"dimk must be 1, 2, or 3, but got {params.dimk}")
    
    # 考虑周期性边界条件下的相邻格点
    for R in Rlist:
        
        # 计算考虑周期性的距离（分数坐标）
        diff = coord2 + R - coord1
        
        # 转换为笛卡尔坐标
        cart_diff = np.dot(diff, poscar_data["lattice"])
        
        # 计算笛卡尔坐标下的距离
        distance = np.linalg.norm(cart_diff)

        if (distance > params.maxdistance) or (distance < params.mindist):
            # 超过最大距离 或者 小于最小距离则跳过，防止出现跃迁到自身的情形
            # print(f"{coord1} to {coord2} distance: {distance} is out of range, skip")
            continue
        
        # 计算耦合强度
        coupling = hopping_strength(distance)
        
        # 存储结果
        coupling_values.append(coupling_sign * coupling)
        distance_values.append(distance)
        R_vectors.append(R)
    
    return coupling_values, R_vectors, distance_values

def calculate_all_couplings(poscar_data, selected_elements, t0=1.0, max_neighbors=1, t0_distance=None, max_distance=10):
    """
    计算两种原子类型之间的所有跃迁强度和向量
    
    参数:
        poscar_data (dict): 通过read_poscar函数读取的结构数据字典
        selected_elements (list): 需要计算耦合的元素列表，如["Mn", "O"]
        t0 (float): 基准跃迁强度，默认为1.0
        max_neighbors (int): 考虑的最大邻居格点数，默认为1
        t0_distances (dict): 原子对的参考距离字典，格式为{(元素1, 元素2): 距离}
        
    返回:
        list: 包含所有耦合信息的列表，每个元素为一个字典，包含:
            - atom1_index: 第一个原子的索引
            - atom2_index: 第二个原子的索引
            - element1: 第一个原子的元素类型
            - element2: 第二个原子的元素类型
            - coupling_values: 耦合强度数组
            - R_vectors: 格矢量数组
    """
    # 获取指定元素类型的原子索引
    all_atom_indices = [i for i, element in enumerate(poscar_data["atom_symbols"]) if element in selected_elements]
    # print(atom_indices1, atom_indices2)
    
    # 存储所有耦合信息
    all_couplings = []

    n_atoms = len(all_atom_indices)
    
    # 计算所有可能的原子对之间的耦合
    for i in range(n_atoms):
        for j in range(i, n_atoms):
            atom1_index = all_atom_indices[i]
            atom2_index = all_atom_indices[j]
            # 跳过相同原子的情况
            # 这里有很大问题，如果是比较远的跃迁，是可以跃迁到自身的
            # if atom1_index == atom2_index:
            #     continue
                
            # 计算耦合强度和格矢量
            coupling_values, R_vectors, distance_values = get_coupling_strength(
                atom1_index, 
                atom2_index, 
                poscar_data
            )
            
            # 存储结果
            coupling_info = {
                "atom1_index": atom1_index,
                "atom2_index": atom2_index,
                "elements": selected_elements,
                "coupling_values": coupling_values,
                "distance_values": distance_values,
                "R_vectors": R_vectors
            }
            
            all_couplings.append(coupling_info)
    
    return all_couplings

def save_hop_to_pythtb(couplings, filename=None):
    """
    将耦合信息转换为pythtb的格式，可以另外保存为文件
    """
    s=""
    for hop in couplings:
        for t, R, d in zip(hop['coupling_values'], hop['R_vectors'], hop['distance_values']):
            s+=f"mymodel.set_hop({t}, {hop['atom1_index']}, {hop['atom2_index']}, [{R[0]}, {R[1]}, {R[2]}]) # distance: {d}\n"
    if filename:
        with open(filename, 'w') as f:
            f.write(s)
    return s

def remove_duplicate_hoppings(couplings):
    """
    消除重复的跃迁对，如果两个跃迁的原子对和R矢量相同或仅相差负号，则认为是重复的跃迁
    
    参数:
        couplings (list): 跃迁信息列表
        
    返回:
        list: 去除重复后的跃迁信息列表
    """
    unique_hoppings = []
    seen_hoppings = set()
    
    for hop in couplings:
        for t, R, d in zip(hop['coupling_values'], hop['R_vectors'], hop['distance_values']):
            # 创建标准化的跃迁标识符
            # 确保原子索引较小的在前面
            if hop['atom1_index'] <= hop['atom2_index']:
                atom_pair = (hop['atom1_index'], hop['atom2_index'])
                R_vec = tuple(R)
            else:
                atom_pair = (hop['atom2_index'], hop['atom1_index'])
                R_vec = tuple(-r for r in R)
                
            hopping_id = (atom_pair, R_vec)
            neg_hopping_id = (atom_pair, tuple(-r for r in R_vec))
            
            # 如果这个跃迁或其反向跃迁都没见过，则添加
            if hopping_id not in seen_hoppings and neg_hopping_id not in seen_hoppings:
                seen_hoppings.add(hopping_id)
                new_hop = {
                    "atom1_index": hop['atom1_index'],
                    "atom2_index": hop['atom2_index'],
                    "elements": hop['elements'],
                    "coupling_values": [t],
                    "distance_values": [d],
                    "R_vectors": [list(R)]
                }
                unique_hoppings.append(new_hop)
    return unique_hoppings

def create_pythtb_model(params=None):
    """
    根据POSCAR文件创建pythtb模型
    
    参数:
        poscar_filename (str): POSCAR文件路径
        params (Parameters): 参数实例，如果为None则使用全局params

    返回:
        pythtb.tb_model: 创建的紧束缚模型
    """
    if params is None:
        params = globals()['params']
        
    try:
        from pythtb import tb_model
    except ImportError:
        raise ImportError("请安装pythtb库: pip install pythtb（建议新建虚拟环境）")
    # print(nspin)
    
    # 读取POSCAR文件
    poscar_data = read_poscar(params.poscar, selected_elements=params.use_elements)
    # print("total_atoms:", len(poscar_data["coords"]))
    
    # 提取晶格向量
    lattice = poscar_data['lattice']/params.a0
    
    # 提取原子坐标（转换为分数坐标）
    coords = poscar_data['coordinates']
    
    # 创建紧束缚模型
    # 3维模型，1个轨道/原子，晶格向量，原子坐标
    model = tb_model(params.dimk, params.dimr, lattice, coords, nspin=params.nspin)
    
    # 计算所有指定元素对之间的耦合
    all_couplings = []
    all_couplings = calculate_all_couplings(
        poscar_data, 
        params.use_elements
    )
    # all_couplings.extend(couplings)
    # print(all_couplings)
    all_couplings = remove_duplicate_hoppings(all_couplings)
    # print(all_couplings)
    
    # 初始化在位能
    for ind in range(len(params.onsite_energy)):
        model.set_onsite(0, ind)
    # 设置磁性
    maglist = params.get_maglist()
    for ind in range(len(maglist)):
        model.set_onsite(maglist[ind]*params.sigma_z, ind, "add")
    # 加上在位能
    for ind in range(len(params.onsite_energy)):
        model.set_onsite(params.onsite_energy[ind]*params.sigma_z, ind, "add")
    
    # 设置跃迁参数
    for hop in all_couplings:
        for t, R, d in zip(hop['coupling_values'], hop['R_vectors'], hop['distance_values']):
            model.set_hop(t, hop['atom1_index'], hop['atom2_index'], [int(R[0]), int(R[1]), int(R[2])], allow_conjugate_pair=True)
            if params.is_print_tb_model_hop:
                print(f"model.set_hop({t}, {hop['atom1_index']}, {hop['atom2_index']}, [{R[0]}, {R[1]}, {R[2]}]) # distance: {d}")

    if params.is_print_tb_model:
        model.display()
    
    return model

def adjust_degenerate_bands(evals, evecs, model, energy_threshold=1e-3):
    """
    修改简并能带的自旋投影。如果找到简并的能带且自旋投影相反，则将其投影设为0。
    
    参数:
        evals (numpy.ndarray): 能带本征值，形状为 (n_bands, n_kpoints)
        evecs (numpy.ndarray): 能带本征矢量，形状为 (n_bands, n_kpoints, n_orbitals, 2)
        model (pythtb.tb_model): 紧束缚模型
        energy_threshold (float): 判断能带简并的能量阈值
        
    返回:
        numpy.ndarray: 修改后的本征矢量
    """
    n_bands, n_kpoints = evals.shape
    modified_evecs = deepcopy(evecs)
    
    # 遍历每个k点
    for k in range(n_kpoints):
        # 遍历每个能带
        for band1 in range(n_bands):
            for band2 in range(band1 + 1, n_bands):
                # 检查能带是否简并
                if abs(evals[band1, k] - evals[band2, k]) < energy_threshold:
                    # 计算两个能带的自旋投影
                    spin_up1, spin_down1 = 0.0, 0.0
                    spin_up2, spin_down2 = 0.0, 0.0
                    
                    for orb in range(model._norb):
                        spin_up1 += abs(evecs[band1, k, orb, 0])**2
                        spin_down1 += abs(evecs[band1, k, orb, 1])**2
                        spin_up2 += abs(evecs[band2, k, orb, 0])**2
                        spin_down2 += abs(evecs[band2, k, orb, 1])**2
                    
                    total1 = spin_up1 + spin_down1
                    total2 = spin_up2 + spin_down2
                    
                    if total1 > 1e-10 and total2 > 1e-10:
                        spin_up1 /= total1
                        spin_down1 /= total1
                        spin_up2 /= total2
                        spin_down2 /= total2
                        
                        # 检查自旋投影是否相反
                        if (abs(spin_up1 - spin_down2) < 0.001 and 
                            abs(spin_down1 - spin_up2) < 0.001):
                            # 将两个能带的自旋投影都设为0
                            for orb in range(model._norb):
                                modified_evecs[band1, k, orb] = 0
                                modified_evecs[band2, k, orb] = 0
    
    return evals, modified_evecs

def check_flat_bands(evals, threshold=0.05, min_points=10):
    """
    检查能带是否存在平带
    
    参数:
        evals (numpy.ndarray): 能带本征值，形状为 (n_bands, n_kpoints)
        threshold (float): 判断平带的能量波动阈值，默认为0.05 eV
        min_points (int): 判断平带的最小连续k点数，默认为10
        
    返回:
        list: 包含平带信息的字典列表，每个字典包含以下键:
            - band_index: 平带所在的能带索引
            - avg_energy: 平带区域的平均能量
            - std_energy: 平带区域的能量标准差
            - k_range: 平带区域的k点范围 [start_k, end_k]
    """
    flat_bands = []
    n_bands, n_kpoints = evals.shape
    
    # 遍历每个能带
    for band in range(n_bands):
        # 寻找连续的平带区域
        start_k = 0
        while start_k < n_kpoints:
            # 找到一段可能的平带区域
            end_k = start_k + 1
            while end_k < n_kpoints:
                # 计算这段区域内能量的标准差
                band_segment = evals[band, start_k:end_k+1]
                std_energy = np.std(band_segment)
                
                # 如果标准差超过阈值，说明不是平带
                if std_energy > threshold:
                    break
                end_k += 1
            
            # 如果找到足够长的平带区域
            if end_k - start_k >= min_points:
                band_segment = evals[band, start_k:end_k]
                flat_bands.append({
                    'band_index': band,
                    'avg_energy': np.mean(band_segment),
                    'std_energy': np.std(band_segment),
                    'k_range': [start_k, end_k]
                })
            
            # 从当前平带区域结束位置继续搜索
            start_k = end_k
            
    return flat_bands



def calculate_band_structure(model, params=None):
    """
    计算能带结构并绘图
    
    参数:
        model (pythtb.tb_model): 紧束缚模型
        params (Parameters): 参数实例，如果为None则使用全局params
    """
    if params is None:
        params = globals()['params']
        
    # 计算能带
    (k_vec, k_dist, k_node) = model.k_path(params.kpath, params.num_k_points)
    evals = model.solve_all(k_vec)
    
    # 调整简并能带
    if params.is_black_degenerate_bands:
        evecs = model.solve_all(k_vec, eig_vectors=True)
        evecs = adjust_degenerate_bands(evals, evecs, model, params.energy_threshold)
    
    # 绘图
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # 绘制能带
    for band in range(evals.shape[0]):
        ax.plot(k_dist, evals[band], 'b-', linewidth=1.0)
    
    # 设置x轴刻度
    ax.set_xticks(k_node)
    ax.set_xticklabels(params.klabel)
    
    # 添加垂直线
    for node in k_node:
        ax.axvline(x=node, linewidth=0.5, color='k')
    
    # 设置标签
    ax.set_xlabel('k-path')
    ax.set_ylabel('Energy (eV)')
    
    # 保存图片
    plt.savefig(f"{params.output_filename}.{params.output_format}")
    plt.close()

def run_band_calculation(poscar_filename):
    """
    运行能带计算
    
    参数:
        poscar_filename (str): POSCAR文件路径
    """
    model = create_pythtb_model(poscar_filename=poscar_filename)
    calculate_band_structure(model=model)

if __name__ == "__main__":
    run_band_calculation(poscar_filename=os.path.join(params.savedir, params.output_filename)) 