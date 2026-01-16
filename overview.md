要件定義書：FARFAR2向け “軽量RNAspider” entanglement boolean filter

1. 背景と目的

1.1 背景
	•	FARFAR2（fragment assembly）では局所的な遷移を繰り返しながらRNA 3D構造を探索するが、探索途中で**不要な絡み（entanglement / knot-like geometry）**が形成されると、以降の探索が非効率・不安定化する可能性がある。
	•	RNAspiderはRNA 3D構造の絡み（entanglement）を検出するための幾何アルゴリズムを提供するが、Rosetta内部に同等機能をそのまま統合するのはコストが高い。

1.2 目的
	•	FARFAR2の遷移判定（accept/reject、またはスコア）に利用可能な、高速な boolean entanglement filter を実装し、探索空間を「絡みの起きにくい領域」に誘導する。

1.3 成果物（Deliverables）
	•	スタンドアロン版：PDB（または座標列）＋二次構造定義 → boolean（entangled / not）を返すフィルタ
	•	Rosetta統合版（段階導入）：
	1.	後処理評価器（trajectory解析）
	2.	soft penaltyとしてのスコア項、またはFilter
	3.	hard reject（必要に応じて）

⸻

2. スコープ

2.1 スコープ内（In scope）
	•	FARFAR2が想定する pseudoknot-free 二次構造制約下での、局所遷移中の絡み抑制。
	•	boolean判定（最初の交差検出で “reject/penalty”）を主目的とする。
	•	closed element（主に loop）に対応する「簡易面」を定義し、backbone線分の貫通を検出する。

2.2 スコープ外（Out of scope）
	•	RNAspiderの完全互換（L/S/D分類、puncture回数・部位の厳密同定）
	•	pseudoknot order / layer再構成（higher-order entanglement探索）
	•	トポロジカルに厳密なノット不変量による保証（本件はヒューリスティクス）

⸻

3. 前提・設計方針

3.1 前提
	•	FARFAR2の二次構造制約（目標二次構造）は信頼でき、探索中も主にその枠内で遷移する。
	•	目的は“探索中の不要な絡みを抑制”であり、偽陽性は許容（探索が多少厳しくなる）し、偽陰性は可能な限り抑える。

3.2 設計方針
	•	Rosetta依存を最小化するため、コアロジックは 入力を単純なデータ構造に限定する。
	•	高速化は「交差判定式」よりも「候補削減（早期終了・近傍化・キャッシュ）」を重視する。
	•	デバッグ性を最優先し、判定根拠（どのループ・どの線分が原因か）を追えること。

⸻

4. 用語定義
	•	Backbone polyline：選択した原子（例：C4′）座標列から作る折れ線（線分集合）。
	•	Closed element：二次構造上、canonical base pairにより閉じられた loop（hairpin/internal/multi等）。
	•	Surface（簡易面）：closed elementに対して定義する面近似。貫通判定の基準。
	•	Entangled（本要件）：backboneの線分が、いずれかの closed element surface を貫通すると判定される状態。

⸻

5. 入出力仕様

5.1 入力（最低限）
	•	座標：
	•	N残基に対して、代表原子座標列（例：C4′、代替可：C1′等）
	•	表現例：vector<Vec3> coords（残基インデックス対応）
	•	二次構造（pseudoknot-free）：
	•	dot-bracket または base-pair list（i,j）の集合
	•	これから closed elements（loop境界）を列挙可能であること
	•	オプション：
	•	今回動いた残基範囲（差分更新のための hint）
	•	判定対象ループ種（hairpin/internalのみ等）

5.2 出力
	•	必須：
	•	bool entangled
	•	推奨（デバッグ用メタ情報、オン/オフ可能）：
	•	loop_id（最初に引っかかった closed element）
	•	segment_id（貫通に寄与した backbone 線分）
	•	近似交点（任意）
	•	判定モード（plane/mesh等）

⸻

6. 機能要件（Functional Requirements）

FR-1：closed elements の構築
	•	二次構造入力（base pairs）から closed elements（loop集合）を構築できること。
	•	最小実装では、以下を対象とする：
	•	hairpin loop
	•	internal loop
	•	multi-loop は初期バージョンでは「任意対応」とし、将来拡張枠に置く。

FR-2：surface（簡易面）生成
	•	closed element に対して、以下いずれかの面近似を生成できること：
	•	Mode A（推奨・最初に実装）：best-fit plane + 2D polygon 内判定
	•	Mode B（オプション）：粗い三角形分割（少数三角形）
	•	surface生成が不安定な場合（点が退化、ほぼ一直線等）には、安全側のフォールバックを定義する：
	•	例：当該ループ判定をスキップ、または保守的に “entangled扱い” など（運用で選択可能）

FR-3：貫通（puncture）判定（boolean）
	•	backboneの線分集合と surface の貫通を判定できること。
	•	boolean filterとして、以下を満たす：
	•	1回でも貫通が見つかったら即座に entangled=true を返す（早期終了）
	•	ループ単位でも早期終了する（当該ループで見つかれば、そのループの残りの判定は不要）
	•	全体として entangled が確定したら終了する

FR-4：候補削減（最低限）
	•	ループごとに bounding volume（AABB等）を保持し、backbone線分のAABBと交差しない場合は判定をスキップできること（broad-phase）。
	•	最初のバージョンは簡易AABBで良い（BVH等は将来拡張）。

FR-5：Rosetta統合（段階導入）
	•	Step 1：後処理（生成構造群への適用）
	•	Step 2：soft penaltyとして統合（スコア項）
	•	Step 3：hard reject（遷移棄却）
	•	実装は「コアロジック＝Rosetta非依存」「アダプタ層＝Rosetta依存」に分割する。

⸻

7. 非機能要件（Non-Functional Requirements）

NFR-1：性能
	•	1遷移あたりの判定時間は、FARFAR2のムーブ頻度に耐えること。
	•	性能目標（暫定）：
	•	cheap gate（AABB）で大半を落とす
	•	expensive判定は少数候補に限定
	•	将来、必要に応じて差分更新・近傍探索を追加できる設計であること。

NFR-2：安定性・ロバスト性
	•	退化ケース（ループ点が少ない、ほぼ共線、座標が欠損）でもクラッシュしない。
	•	数値誤差に対する許容幅（epsilon）をパラメータ化し、再現性を確保する。

NFR-3：デバッグ容易性
	•	判定根拠を出力できる（ログ、あるいは戻り値のメタ情報）。
	•	再現用に、該当構造をダンプするフック（任意）を用意可能にする。

NFR-4：保守性
	•	Rosetta側変更に引きずられないよう、入力インタフェースを最小化する。
	•	コアロジックは単体テスト可能であること。

⸻

8. インタフェース設計（モジュール構造）

8.1 コアライブラリ（Rosetta非依存）
	•	EntanglementFilterCore
	•	build_loops(base_pairs, n_res) -> vector<Loop>
	•	precompute_surfaces(coords, loops) -> vector<Surface>
	•	is_entangled(coords, loops/surfaces, options) -> Result

8.2 Rosettaアダプタ（Rosetta依存）
	•	Pose から座標列（代表原子）を抽出
	•	FARFAR2が保持する二次構造制約（あるいは目標二次構造）から base-pair list を取得
	•	遷移判定の統合点へ boolean / penalty を返す

⸻

9. パラメータ（初期値は暫定）
	•	代表原子：C4′（推奨）、代替：C1′
	•	surface mode：A（plane + polygon）をデフォルト
	•	epsilon：
	•	平面跨ぎ判定閾値
	•	点内判定の境界許容
	•	対象ループ種：hairpin/internal（初期）
	•	判定頻度：毎ムーブ（初期）、将来 Nステップごとに変更可能

⸻

10. テスト計画（受け入れ基準）

10.1 単体テスト
	•	既知の簡単構造（絡まない）で entangled=false
	•	人工例（明確に貫通する）で entangled=true
	•	退化ケースでクラッシュしない

10.2 統合テスト（後処理）
	•	FARFAR2で生成した構造集合に適用し、
	•	entangled判定率
	•	どのループが原因になりやすいか
	•	判定コスト分布
をレポートできること

10.3 Rosetta統合テスト
	•	soft penalty統合で動作する（クラッシュしない、ログが出る）
	•	hard rejectを有効化しても探索が停止しない（過度な偽陽性がない）
※過度な場合は multi-loop除外や判定頻度低減で調整可能であること

⸻

11. 開発計画（マイルストーン）
	1.	M1：スタンドアロン core（plane近似＋早期終了＋簡易AABB）
	2.	M2：後処理評価（FARFAR2出力に適用、統計・デバッグ）
	3.	M3：Rosetta soft penalty 統合（最小差分）
	4.	M4：性能改善（差分更新、候補削減強化）
	5.	M5：必要なら hard reject、multi-loop対応

⸻

12. 未決事項（最初に決めると良い項目）
	•	multi-loop の扱い：初期は除外か、分割して扱うか
	•	surface生成が不安定なループのフォールバック方針：skipかrejectか
	•	代表原子：C4′固定か、状況により切替可能にするか
	•	統合点：acceptance直前か、ムーブ直後か、minimization前後か

⸻

必要なら、この要件定義書をベースに、次の一手として「M1（スタンドアロン core）」のAPIシグネチャとデータ構造（Loop/Surface/Result）、および「plane近似モード」の判定手順を、実装に落ちる粒度まで具体化します。