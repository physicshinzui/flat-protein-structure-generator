import os
import numpy as np
import pymol
from pymol import cmd
from tqdm import tqdm

pymol.finish_launching(['pymol', '-c'])

def build_flat_peptide_and_save(peptide_seq, data_dict):
    """
    PyMOLのfabコマンドを使って平坦なペプチド構造を生成し、
    原子座標などをNumPyの辞書に保存する。

    Parameters
    ----------
    peptide_seq : str
        1文字表記のアミノ酸配列 (例: "ACDEFGHIK")
    data_dict : dict
        生成したデータを格納する辞書
    """

    # PyMOLの状態をリセット（既存オブジェクトを消去）
    cmd.reinitialize()

    # 1) fabコマンドで"flat"ペプチドを生成
    cmd.fab(peptide_seq, "peptide", ss=4)

    # 2) PyMOLから原子モデルを取得
    atoms = cmd.get_model("peptide")

    coords_list = []
    atom_names_list = []
    residue_ids_list = []
    residue_names_list = []

    # 3) 取得したモデルの各原子情報をリスト化
    for a in atoms.atom:
        coords_list.append([a.coord[0], a.coord[1], a.coord[2]])  # (x,y,z)
        atom_names_list.append(a.name)    # 例: 'N', 'CA', 'C', 'O' など
        residue_ids_list.append(a.resi_number)  # 数値としての残基番号 (1,2,3,...)
        residue_names_list.append(a.resn)       # 3文字表記: 'ALA', 'VAL', ...

    # 4) NumPy配列に変換
    coordinates_array = np.array(coords_list, dtype=np.float32)  # shape=(N,3)
    atom_names_array = np.array(atom_names_list, dtype='S6')  # shape=(N,)
    residue_ids_array = np.array(residue_ids_list, dtype=np.int32)  # shape=(N,)
    residue_names_array = np.array(residue_names_list, dtype='S6')  # shape=(N,)

    # 5) データを辞書に格納
    data_dict[peptide_seq] = {
        "coordinates": coordinates_array,
        "atom_names": atom_names_array,
        "residue_ids": residue_ids_array,
        "residue_names": residue_names_array
    }


if __name__ == "__main__":
    nskips = 10
    output_npz = f"flat_structures_{nskips}.npz"
    if os.path.exists(output_npz):
        os.remove(output_npz) 

    peptide_list = np.load('../all_training_seqs.npy')[::nskips]
    #peptide_list = np.load('../v2/all_training_seqs_skip10000.npy')

    data_dict = {}
    print(peptide_list)
    for seq in tqdm(peptide_list):
        build_flat_peptide_and_save(seq, data_dict)

    # NumPyの辞書形式に変換し、順序情報も追加
    save_dict = {}
    for peptide, data in data_dict.items():
        for key, value in data.items():
           save_dict[f"{peptide}_{key}"] = value
    
    # 順序情報を保存（後で復元可能）
    save_dict["sequence_order"] = np.array(peptide_list, dtype="S")

    # .npzファイルに保存
    np.savez(output_npz, **save_dict)

    print(f"生成完了: {output_npz}")

