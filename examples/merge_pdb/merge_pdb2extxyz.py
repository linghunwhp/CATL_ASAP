#!/usr/bin/env python3
# merge_pdb2extxyz.py
import argparse, glob, os, ase.io
from ase.io import extxyz

def merge_pdb_to_extxyz(pattern, outfile, multiframe=False):
    """
    pattern : 通配符，如 '*.pdb'
    outfile : 输出 extxyz 文件名
    multiframe: False -> 场景 A（单帧）
                True  -> 场景 B（多帧轨迹）
    """
    pdb_files = sorted(glob.glob(pattern))
    if not pdb_files:
        raise RuntimeError(f'未找到匹配 {pattern} 的 PDB 文件')

    all_images = []
    for pdb in pdb_files:
        # ASE 读 PDB 会自动识别 CRYST1 记录
        atoms = ase.io.read(pdb)
        # 把 CONECT 转成 extxyz 格式
        conect = atoms.info.get('conect', {})
        conn_lines = [" ".join(str(x) for x in [i]+partners)
                      for i, partners in conect.items()]
        atoms.info['connectivity'] = "\n".join(conn_lines)
        all_images.append(atoms)

    if multiframe:
        # 场景 B：直接写多帧
        extxyz.write_extxyz(outfile, all_images, write_info=True)
    else:
        # 场景 A：合并成单帧
        # 1. 取第一个文件的 Cell 作为统一 Cell
        cell = all_images[0].get_cell()
        # 2. 把所有坐标拼到一起
        symbols = []
        positions = []
        info = {'connectivity': []}
        offset = 0   # 原子全局序号偏移
        for atoms in all_images:
            symbols.extend(atoms.get_chemical_symbols())
            positions.extend(atoms.get_positions())
            # 调整 connectivity 全局序号
            for i, partners in atoms.info.get('conect', {}).items():
                info['connectivity'].append(
                    " ".join(str(x+offset) for x in [i]+partners))
            offset += len(atoms)
        # 3. 新建一个 Atoms
        big_atoms = ase.Atoms(symbols=symbols, positions=positions, cell=cell, pbc=True)
        big_atoms.info['connectivity'] = "\n".join(info['connectivity'])
        extxyz.write_extxyz(outfile, big_atoms, write_info=True)

    print(f'已写入 {outfile} ，共 {len(all_images)} 个 PDB 来源。')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='把多个 PDB 合并成单个 extxyz')
    parser.add_argument('-i', '--input', default='*.pdb',
                        help='PDB 通配符 (默认: *.pdb)')
    parser.add_argument('-o', '--output', default='merged.xyz',
                        help='输出 extxyz 文件名 (默认: merged.xyz)')
    parser.add_argument('-m', '--multiframe', action='store_true',
                        help='生成多帧轨迹，而非合并成单帧')
    args = parser.parse_args()

    merge_pdb_to_extxyz(args.input, args.output, multiframe=args.multiframe)