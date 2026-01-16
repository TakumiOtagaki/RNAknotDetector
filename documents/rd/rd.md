要件定義書：RNA Loop–Backbone Entanglement Penalty（C4′, plane model）

0. 目的

与えられた RNA 構造（C4′座標）と二次構造（loops）に対し、各ループが張る近似面を backbone が貫通する回数（またはイベント数）を評価し、MC遷移での受理確率を 増分 \Delta K に基づき抑制するためのスカラー量を返す。

⸻

1. 入力・出力

1.1 入力
 - loops : [ (looptype, [(i, j), ...]) ]
 - looptype ∈ {hairpin, interior, stacking, multiloop}
 - [(i,j), ...] はそのループ（構造要素）を閉じる base pairs（インデックスは residue index）
 - hairpin: 1組
 - interior/stacking: 2組
 - multiloop: 3組以上
 - r_{i,k} = (x_{i,k}, y_{i,k}, z_{i,k})
 - residue i の原子 k の3D座標
 - 本仕様では k = "C4'" のみを使用する（欠損時の扱いはNFRで規定）

1.2 出力
 - K：entanglement イベント数（整数）または重み付きスコア（実装方針により選択）
 - details（任意、デバッグ用）：
 - 最初に貫通した (loop_id, segment_id) 等

⸻

2. 定義

2.1 Backbone polyline
 - 点列：p_i = r_{i,"C4'"}
 - 線分：s_i = (p_i, p_{i+1}) for i = 1..n-1（存在する場合のみ）

2.2 ループの“所属残基集合” R_l

ループ l に対し、「そのループを構成する残基集合」 R_l を定義する。
これは 貫通判定から除外する残基集合である。
 - hairpin (closing pair = (i,j)):
 - R_l = \{i, i+1, ..., j\}
 - interior/stacking (closing pairs = (i1,j1), (i2,j2)):
 - 2本の境界鎖を持つ。典型的には
 - chain A: [i1, ..., i2]
 - chain B: [j2, ..., j1]（向きは実装上は集合なので無視）
 - R_l = [i1..i2] \cup [j2..j1]
 - multiloop (closing pairs = (a1,b1), (a2,b2), …, (am,bm)):
 - 各 branch の端点ペアを与えられていると解釈し、全 branch を包含する “外周” は曖昧になりやすい。
 - 初期実装では、branchごとに独立 surface を作るため、各 branch t について
 - R_{l,t} = [a_t..a_{t+1}] \cup [b_{t+1}..b_t] などの局所領域
 - 実装簡略化のため、最初は multiloop を interior の集合（複数面）として扱う（FRで明示）。

注：上の区間の並びは、[(i,j), ...] の並び順（5’→3’）と整合することを仮定します。整合が取れない場合は、追加の “loop boundary extraction” ロジックが必要です（本書では要件外）。

2.3 貫通判定から除外する線分（スキップ規則）

ループ l の判定で backbone 線分 s_u = (p_u, p_{u+1}) を評価する際、次を満たすならスキップ：
 - u ∈ R_l または u+1 ∈ R_l
（＝線分の端点がループ構成残基に触れていれば除外）

これにより boundary/closing 近傍の数値誤差由来の偽陽性を抑える。

⸻

3. ループ surface の近似（plane model）

各ループ（または multiloop の各 branch）に対し、定数個 C の平面を作る。最初は C=1。

3.1 平面の定義

平面は (c, n_hat) で表す。
 - c：平面上の点（通常は点群重心）
 - n_hat：単位法線ベクトル

3.2 best-fit plane の計算

ループ境界点集合 B_l = { p_t | t ∈ boundary_indices(l) } から計算：
	1.	c = mean(B_l)
	2.	共分散 Σ = Σ_t (p_t - c)(p_t - c)^T
	3.	n_hat は Σ の最小固有値に対応する固有ベクトル（点群の最小分散方向）

退化（|B_l| < 3、ほぼ共線など）の場合：
 - surface_valid=false とし、当該ループは評価から除外（または保守的にペナルティ）
初期は 除外推奨（探索を殺しにくい）。

3.3 平面上での“領域内判定”

ループ境界点 B_l を平面へ射影し、2D多角形として扱う。
 - 平面基底を作る：
 - e1：n_hat に直交する任意単位ベクトル
 - e2 = n_hat × e1
 - 射影：
 - 各点 q_t = ( dot(p_t - c, e1), dot(p_t - c, e2) )（2D）
 - P_l = [q_t] を多角形（頂点列）とみなし、点 q が内部かを判定する。

boolean filter用途で高速化したい場合は、P_l を
 - 三角形分割（扇形）して定数枚にする、または
 - 凸包近似にする
などで内判定を定数化してよい。

⸻

4. 線分が surface を“貫通”する判定

ループ l の平面 π=(c,n_hat) と backbone 線分 s=(a,b) について：

4.1 平面跨ぎ（交差候補）

d_a = \langle a-c, n\_hat \rangle,\quad d_b = \langle b-c, n\_hat \rangle
 - d_a * d_b > 0：同じ側 → 交差しない（スキップ）
 - |d_a| < eps または |d_b| < eps：境界近傍 → 初期実装ではスキップ（偽陽性抑制）
 - それ以外：交点を計算

4.2 交点

t = \frac{d_a}{d_a - d_b}\quad (0<t<1)
x = a + t(b-a)

4.3 領域内判定

交点 x を 2D に射影：q = ( dot(x-c,e1), dot(x-c,e2) )
q ∈ interior(P_l) なら “puncture event” とカウント。

⸻

5. 全体アルゴリズム（K の計算）

5.1 フル評価（1構造に対する K）

Input: loops, r_{i,k}, n
Output: K

1. p_i <- r_{i,"C4'"} for i=1..n (skip residues with missing C4')
2. segments S <- { (i, p_i, p_{i+1}) | i=1..n-1 and both points exist }

3. K <- 0
4. for each loop l in loops:
     (a) determine boundary indices and residue set R_l
     (b) build plane(s) Π_l = {π_1,...,π_C} using boundary points
     (c) build 2D polygon(s) P_l for each plane (projection)

     for each segment s_i in S:
        if i in R_l or (i+1) in R_l: continue  // skip rule
        // (optional) AABB/broad-phase check; see Section 6
        for each plane π in Π_l:
            if SegmentIntersectsLoopSurface(s_i, π, P_l):
                K <- K + 1
                // early-exit policy (optional):
                // if boolean-only: return K>0
                // if counting: can break inner loops or continue, depending on definition
5. return K

5.2 増分評価（MCに使う ΔK：推奨）

MC遷移で変更された残基集合 M（フラグメント区間など）が分かるなら：
 - 評価対象の線分集合を
S_delta = { s_i | i ∈ [min(M)-1, max(M)] } に限定
 - ループ側も、AABBや残基距離で候補ループだけに絞る

そして K_new_on_delta と K_old_on_delta の差から \Delta K を得る。

初期実装では、まずフル評価で正しさ確認 → 次に差分へ。

⸻

6. AABB とは何か（あなたの質問への回答）

**AABB = Axis-Aligned Bounding Box（軸平行バウンディングボックス）**です。
3D空間で「ある集合（点群、線分、面など）を全部含む、xyz軸に平行な直方体」を指します。
 - ある点集合 X の AABB は
[\min x, \max x] \times [\min y, \max y] \times [\min z, \max z]
 - 2つのAABBが交差しないなら、その中身（線分や面）同士も交差し得ないので、詳細な交差判定をスキップできます。

この “先に粗く弾く” ステップを **broad-phase（候補削減）**と呼び、計算量の定数因子を大きく下げます。

AABB を本アルゴリズムでどう使うか（推奨）
 - 各ループ surface（境界点群）について AABB を事前計算
 - 各 backbone 線分についても AABB は即時計算できる（端点2つの min/max）
 - 交差判定前に
 - AABB(loop) ∩ AABB(segment) == ∅ なら continue

これだけでも高速化になります。

⸻

7. 非機能要件（NFR）
 - 欠損 C4′ がある残基は、当該点を使う線分を生成しない
 - 数値安定性：
 - eps（例：1e-3〜1e-2 Å相当）を導入し、平面境界付近の判定をスキップ
 - デバッグ：
 - どのループ・どの線分でヒットしたかログ可能

⸻

8. 未決事項（実装前に決めるとよい）
	1.	loops の (i,j) 列の順序保証（5’→3’で整っているか）
	2.	multiloop の分割規則（branchごとに interior と同様に作る、でよいか）
	3.	K の定義：
 - 「ループ×線分」で1カウント
 - 「ループごとに最大1カウント」
 - など（あなたは soft penalty 用なので、まずは “新規ヒット数” で良い）

⸻

次のステップ（実装に入るための最短）

あなたが実装に移るなら、まずは Section 5.1 のフル評価を Python か C++ で作り、少数PDBで挙動確認。その後、FARFAR2統合に合わせて Section 5.2 の差分評価＋AABBを入れる、が一番堅いです。

必要なら、上記仕様に沿って「関数分割（インタフェース）」まで提案します（例：build_loop_boundary_indices, fit_plane, project_polygon, segment_plane_intersection, point_in_polygon_2d, aabb_intersect など）。