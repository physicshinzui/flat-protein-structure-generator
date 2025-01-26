import numpy as np

def npz_to_pdb(npz_path, peptide, pdb_output_path, chain_id='A'):

    data = np.load(npz_path, allow_pickle=True)

    coord_key = f"{peptide}_coordinates"
    atom_key  = f"{peptide}_atom_names"
    resi_key  = f"{peptide}_residue_ids"
    resn_key  = f"{peptide}_residue_names"

    coordinates = data[coord_key]
    atom_names  = data[atom_key]
    residue_ids = data[resi_key]
    residue_names = data[resn_key]

    # バイト列を文字列に (必要なら)
    atom_names = [a.decode() if isinstance(a, bytes) else str(a) for a in atom_names]
    residue_names = [r.decode() if isinstance(r, bytes) else str(r) for r in residue_names]

    # ▼ PDBフォーマット (ATOM行) の作り方
    # 各カラム位置を厳密に定義
    #
    #  Format: 
    #  1-6    6s    : レコード名 "ATOM  "
    #  7-11   5d    : 原子番号 (右詰め)
    #  12     1s    : 空白
    #  13-16  4s    : 原子名 (場所に注意)
    #  17     1s    : 可替位置 (altLoc) 今回は空白
    #  18-20  3s    : 残基名 (右詰め)
    #  21     1s    : 空白
    #  22     1s    : チェーンID
    #  23-26  4d    : 残基番号
    #  27     1s    : インサートコード 今回は空白
    #  28-30  3s    : 空白
    #  31-38  8.3f  : x座標 (右詰め)
    #  39-46  8.3f  : y座標 (右詰め)
    #  47-54  8.3f  : z座標 (右詰め)
    #  55-60  6.2f  : occupancy
    #  61-66  6.2f  : tempFactor
    #  67-76  10s   : 空白
    #  77-78  2s    : 元素記号 (右詰め)
    #  79-80  2s    : 電荷 今回は空白
    
    # Pythonのf文字列で厳密に揃えづらい場合、formatや%を使う方法もある
    line_format = (
        "{record:<6}"       # Columns 1-6:  "ATOM  " or "HETATM"
        "{atom_id:5d}"      # Columns 7-11: 原子番号
        " "                 # Column 12:    空白
        "{atom_name:^4}"    # Columns 13-16:原子名を中央寄せ (一般的なタンパク質表記)
        "{altLoc:1s}"       # Column 17:    可替位置 (今回は空白に)
        "{res_name:>3}"     # Columns 18-20:残基名 (右詰め)
        " "                 # Column 21:    空白
        "{chain_id:1s}"     # Column 22:    チェーンID
        "{res_id:4d}"       # Columns 23-26:残基番号 (右詰め)
        "{iCode:1s}"        # Column 27:    インサートコード (空白)
        "   "               # Columns 28-30: 空白
        "{x:8.3f}"          # Columns 31-38: x座標
        "{y:8.3f}"          # Columns 39-46: y座標
        "{z:8.3f}"          # Columns 47-54: z座標
        "{occ:6.2f}"        # Columns 55-60: occupancy
        "{temp:6.2f}"       # Columns 61-66: tempFactor
        "          "         # Columns 67-76: 空白
        "{element:>2s}"     # Columns 77-78: 元素記号 (右詰め)
        "{charge:>2s}"      # Columns 79-80: 電荷 (今回は空白)
    )

    with open(pdb_output_path, "w") as f:
        for i, (xyz, aname, rname, rid) in enumerate(zip(coordinates, atom_names, residue_names, residue_ids), start=1):
            x, y, z = xyz
            # 原子名が1文字の場合など、Columns 13-16との兼ね合いに注意
            # 例: aname="N", "CA" などで長さが違う
            # 元素記号っぽい推定 (C, N, O, etc.) するなら先頭1-2文字を見るなど
            # ここでは単純に aname の先頭1-2文字を元素と仮定
            element_guess = aname[0].upper() if aname else "C"

            line = line_format.format(
                record="ATOM",
                atom_id=i,
                atom_name=aname[:4],  # 4文字まで
                altLoc=" ",
                res_name=rname[:3],   # 3文字まで
                chain_id=chain_id,
                res_id=int(rid),
                iCode=" ",
                x=x, y=y, z=z,
                occ=1.00,
                temp=0.00,
                element=element_guess.rjust(1),  # 一文字元素なら右詰
                charge="  "
            )
            f.write(line + "\n")
        f.write("END\n")

    print(f"PDBファイルを出力しました: {pdb_output_path}")

if __name__ == "__main__":
    npz_path = "flat_structures.npz"  # NumPyファイル
    peptide = "AECR"                  # 取得するペプチド
    pdb_output_path = f"{peptide}.pdb" # 出力PDBファイル名

    npz_to_pdb(npz_path, peptide, pdb_output_path, chain_id='A')
