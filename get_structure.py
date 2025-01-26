import numpy as np 

def get_structure_from_npz(npz_path, peptide):
    """
    .npzファイルから指定したペプチドの構造データ(座標, 原子名, 残基名など)を読み込む

    Parameters
    ----------
    npz_path : str
        NumPyの`.npz`ファイルのパス
    peptide : str
        取得したいペプチドのアミノ酸配列 (例: "AAA")

    Returns
    -------
    tuple
        (coordinates, atom_names, residue_ids, residue_names)
    """
    data = np.load(npz_path, allow_pickle=True)

    # 取得
    coordinates = data[f"{peptide}_coordinates"]
    atom_names = data[f"{peptide}_atom_names"]
    residue_ids = data[f"{peptide}_residue_ids"]
    residue_names = data[f"{peptide}_residue_names"]

    return coordinates, atom_names, residue_ids, residue_names



# 使用例
npz_path = "flat_structures.npz"
peptide = "AECR"  # 例: AAAのデータを取得

coordinates, atom_names, residue_ids, residue_names = get_structure_from_npz(npz_path, peptide)

print(f"Coordinates shape: {coordinates.shape}")
print(f"Atom Names: {atom_names}")
print(f"Residue IDs: {residue_ids}")
print(f"Residue Names: {residue_names}")


def get_sequence_order(npz_path):
    """ .npzファイルから保存されたペプチドの順番を取得 """
    data = np.load(npz_path, allow_pickle=True)
    return [s.decode() for s in data["sequence_order"]]

# 順序を取得
sequence_order = get_sequence_order("flat_structures.npz")
print(sequence_order)  # ["AAA", "KKK", "PPP", ...] のように順番通り出力

