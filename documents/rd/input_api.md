入力API要件（案）

目的
- コアロジックはC++で実装し、座標配列と二次構造を受け取る。
- PyMOLからの利用はPython経由を前提とし、BiopythonでPDB/mmCIFを処理する。
- FARFAR2統合時も同じC++コアを使える形にする。

対象範囲
- Python I/O層: PDB/mmCIFの読み込み、代表原子座標の抽出、二次構造の受け取り。
- C++コア層: 座標配列と二次構造からエンタングル判定に必要な入力構造体を作る。

基本フロー（PyMOL想定）
1) PythonでPDB/mmCIFを受け取る
2) Biopythonで指定原子の座標を抽出する
3) PythonからC++コアに座標配列と二次構造を渡す

座標入力の方針
- 代表原子はC4'をデフォルトとする。
- atom_namesで任意の原子名を指定可能にする（例: {"C4'"}）。
- 1残基あたり1原子を基本とし、必要なら2原子まで拡張できる設計とする。

座標データ形式（C++側の案）
struct ResidueCoord {
  int res_index;               // 1-based
  std::vector<Vec3> atoms;     // atom_namesの順序に対応
};

二次構造入力
- dot-bracket
  - 例: "(((...)))"
  - pseudoknot-free を前提
- base-pair list
  - 例: [(i, j)] or [(i, j, bp_type)]
  - bp_typeは任意（例: "WC", "GU", "NONCANON"）
  - bpseq由来の変換を許容

API入力（C++側の案）
struct InputSpec {
  std::vector<ResidueCoord> coords;   // 必須
  std::string dot_bracket;            // 任意
  std::vector<BasePair> bp_list;      // 任意
  std::vector<std::string> atom_names; // 省略時は{"C4'"}
};

struct BasePair {
  int i, j;                // 1-based
  std::string bp_type;     // optional
};

動作ルール
- 二次構造は以下の優先順:
  1) bp_list
  2) dot_bracket（bp_listに変換して扱う）
- coordsが空の場合はエラーとする（座標抽出はPython側の責務）。

デバッグ・可視化
- 判定根拠としてloop_id、segment_id、平面情報（法線+一点）などを返せるようにする。
- PyMOL向けに線分・平面ポリゴン・交点を出力する補助関数を用意する。
- 出力形式は後で決める（例: PDB/XYZ/CGO）。

未決事項
- 複数原子（2点）を使う場合のコアロジック側の扱い
- bp_typeをコアロジックで使うかどうか
- TODO: Python -> C++ で bp_type（canonical/non-canonical）を渡す手段と、その受け口APIを定義する
