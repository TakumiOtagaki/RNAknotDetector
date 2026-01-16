入力API要件（案）

目的
- スタンドアロンでPDB/mmCIFから座標を抽出し、二次構造情報と合わせてコアロジックに渡す。
- FARFAR2統合を見据えて、C++で使える薄いAPIを定義する。
- PyMOLで呼び出す運用も想定し、Python経由の入出力パスを用意する。

対象範囲
- 入力: PDB/mmCIFファイル、または既に抽出済みの座標配列
- 二次構造: dot-bracket または base-pair list（bpseq相当）
- 出力: コアロジックへの入力構造体 + デバッグメタ情報

座標入力の考え方
- 代表原子は基本 C4' をデフォルトとする。
- 任意に原子名を指定可能にする（例: C4', C1'）。
- 座標保持は「残基iと原子種kに対応する座標」を扱える形にする。

座標データ形式（内部表現の案）
- 1残基あたり1原子の代表座標を基本とする。
- 必要なら複数原子（最大2つ程度）に拡張できる設計とする。

例: 代表1原子
struct AtomCoord {
  int res_index;     // 1-based
  int atom_index;    // 0..M-1 (指定原子の順序)
  double x, y, z;
};

例: 残基単位
struct ResidueCoord {
  int res_index;     // 1-based
  std::vector<Vec3> atoms; // atom_indexに対応
};

二次構造入力
- dot-bracket 形式
  - 例: "(((...)))"
  - pseudoknot-free を前提
- base-pair list
  - 例: [(i, j)] or [(i, j, bp_type)]
  - bp_type は任意（例: "WC", "GU", "NONCANON"）
  - bpseq 由来の変換を許容

API入力（C++の案）
struct InputSpec {
  std::string pdb_or_cif_path;     // 任意
  std::vector<ResidueCoord> coords; // 直接入力がある場合
  std::string dot_bracket;         // 任意
  std::vector<BasePair> bp_list;   // 任意
  std::vector<std::string> atom_names; // 例: {"C4'"}、省略時はC4'
};

BasePair:
struct BasePair {
  int i, j;                  // 1-based
  std::string bp_type;       // optional
};

動作ルール
- 座標は以下の優先順:
  1) coords が与えられている場合はそれを使う
  2) それ以外は pdb_or_cif_path から抽出
- 二次構造は以下の優先順:
  1) bp_list
  2) dot_bracket
- dot_bracket は bp_list に変換して扱う

PDB/mmCIFの取り扱い
- PyMOLから呼び出す場合は以下の流れを基本とする:
  1) PythonでPDB/mmCIFを受け取る
  2) Biopythonで座標を抽出する
  3) C++のコアロジックに渡す
- C++側は「座標配列+二次構造」を受け取る形に寄せる。

デバッグ・可視化
- 判定根拠として「loop_id」「segment_id」「平面（法線+一点）」などを返せるようにする。
- PyMOL向けに「線分」「平面ポリゴン」「交点」を出力できる補助関数を用意する。
- 出力形式は後で決める（例: PDB/XYZ/CGO）。

未決事項
- 複数原子（2点）を使う場合のコアロジック側の取り扱い
- bp_type をコアロジックで使うかどうか
- mmCIFの具体的な抽出実装（C++パーサ選定）
