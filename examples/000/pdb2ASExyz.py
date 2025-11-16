# pdb2extxyz.py
import os
import ase.io
from ase import Atoms
from ase.io import extxyz


# 2. 用 ASE 读 PDB（会自动识别晶胞）
in_path = 'pdbbank'  # 假设 PDB 文件都在 pdbbank 目录下
out_path = 'extxyz_output'
file_list = [f for f in os.listdir(in_path) if f.endswith('.pdb')]
for file in file_list:
    atoms = ase.io.read(os.path.join(in_path, file))

    # 3. 把 CONECT 信息转成 extxyz 的 connectivity 字段
    #    ASE 的 PDB 解析器会把 CONECT 记录在 atoms.info['conect']
    conect = atoms.info.get('conect', {})
    # 转成 extxyz 需要的格式：每行 “i j k …” 表示 i 原子与 j k … 相连
    conn_lines = []
    for i, partners in conect.items():
        line = " ".join(str(x) for x in [i] + partners)
        conn_lines.append(line)
    atoms.info['connectivity'] = "\n".join(conn_lines)

    # 4. 写出 extxyz
    extxyz.write_extxyz(os.path.join(out_path, file.replace('.pdb', '.xyz')), atoms, write_info=True)

    print(f"已写入 {file.replace('.pdb', '.xyz')}")