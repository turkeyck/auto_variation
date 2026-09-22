# 規格整理 v 1.2.0

> 產生日期：2026-09-17。本規格基於對 `cadabra/`（4 支既有檔案 + 1 支新 lint script）與
> `sympy_layer/`（12 支檔案，約 2500 行）的逐檔閱讀，以及對 `git log`（11 筆 commit）、
> 目前未 commit 的 diff、執行環境（Windows `py -3` + WSL Ubuntu `cadabra2 2.5.14`）的實測驗證所寫成。
> 範圍已由你確認：**涵蓋 A–J 全部模組的完整長期系統**；cadabra2 執行環境已完成版本調查（見下）。

---

## 0）cadabra2 版本調查結果（依你的要求，先於技術規格完成）

調查方式：讀 `github.com/kpeeters/cadabra2` 的 releases 列表、changelog、issue tracker（issue #166、#387）、
以及 2025-07-31 release 之後的 commit 紀錄（最新一筆到 2025-09-24，全為建置/CI 修正，與本專案無關）。

```md
[調查結論]
- 最新穩定版：2.5.14（2025-07-31 發布），也是目前 WSL Ubuntu 已安裝的版本（`dpkg -l` 實測確認）。
- 2025-09-24 之後 master 分支沒有任何新 release；近期 commit 全部是建置系統（boost、appimage、windows installer）修正，
  與索引檢查/dummy 重新命名/substitute 邏輯無關。
- cadabra_utils.py 記載的三個 bug（Accent 邊界看不穿自由指標檢查、rename_dummies 在 InverseMetric/Symmetric
  縮併上報錯、substitute() RHS 重用啞指標名衝突）在官方 issue tracker 中**沒有對應的已回報 issue**，
  也沒有任何 changelog 條目聲稱修復過這三個行為。
- 唯一相關的既有 issue 是 #166（"Ex_comparator::equal_subtree 未正確處理 use_props"，已於 2.2.6 / 2019-04-16
  修復）——但那是 Accent 比對邏輯的另一個、更早、已修復的問題，不是本專案發現的三個 bug。
- 另一個開放中的 issue #387（"rename_dummies 應改用位置排序而非名稱排序"，2026-01-19 由維護者本人開立，
  目前 open）與本專案的三個 bug 也無直接關聯。
```

```md
[結論與建議]
- 沒有任何目前存在的 cadabra2 版本（含最新 2.5.14 與 master 分支）修復這三個 bug。這不是本機安裝過舊，
  是上游尚未觸及這個行為。
- 建議：(1) 技術規格內明訂固定使用 WSL Ubuntu 內的 cadabra2 2.5.14，並把 D 模組的三條 workaround 規則
  視為長期營運規則，不是等待升級後可拿掉的暫時措施；(2) 建議把這三個 bug 各自整理成最小可重現範例，
  向 kpeeters/cadabra2 開 issue——即使不指望短期修復，也能讓未來每次重新驗證有明確的比對基準；
  (3) 每次 cadabra2 有新 release 時，只需要重跑 `cadabra_utils.py` 文件開頭建議的三個最小重現腳本，
  不需要每次重新驗證全部程式碼。
```

此調查結果已併入下方技術規格的「執行環境」與「開發應注意重點」章節。

---

## 技術規格文件

### 假設與前提

- 本專案沒有圖形介面、沒有 Web 服務、沒有資料庫、沒有檔案上傳、沒有即時對話式 LLM 生成內容。
  它是一套以 Python（cadabra2 + SymPy）撰寫的**符號計算研究程式庫**，由研究者（你）與程式碼協作
  agent（Claude Code / Codex）透過命令列腳本操作，輸出是終端機文字報告與（未來）LaTeX 片段。
  因此本規格中「UI 設計／色彩策略／SVG 畫面示意／task model／state model（UI 版）／UI 元件清單／
  UI 事件回報／UI↔API Mapping／上傳預覽／縮放比例／Streaming／通知與背景執行（UI 意義上的）」
  這些章節，改以本領域對應的等價章節取代，並在對應小節註明取代理由，而非留空或硬套。
- 假設執行環境維持現況：SymPy 層在 Windows 用 `py -3`（已確認 sympy 1.14.0 可用）；cadabra 層在
  WSL Ubuntu 用系統套件 `cadabra2 2.5.14`（已透過 `dpkg -l` 確認版本、透過三個既有 bug 的重新驗證
  紀錄確認可用）。這兩個直譯器目前無法在同一個 Python 環境下共用，因此 cadabra 層與 SymPy 層是
  兩個獨立可執行單元，靠 `docs/SPEC.md` 與（新增的）`verification/` 模組互相對齊，不是靠 import 互通。
- 假設「已完成」的定義採你在需求裡訂下的標準：對照哪一條論文式子、用哪一條獨立路徑交叉驗證，
  兩者缺一都不算完成。本規格中「現況」一律標注這兩項；缺任何一項的既有程式碼，即使目前
  `run_checks()` 顯示 PASS，也在此規格中重新分類為「部分完成」或「需重做」。
- 假設你會持續親自審閱重大推導結果（規格中「我對工作方式的要求」第 2–4 點），因此本規格把
  「產出讓你能快速核對的中間結果」當成一等公民需求，而不只是把最終答案印出來。

### 背景

專案目前有 11 筆 commit，最新一筆完成了廣義 Proca 背景理論（對照 arXiv:1603.05806）加入 L5 項，
6/7 項自我檢查通過（已重新實測確認），1 項已知問題誠實標記為未解決。現有程式碼呈現一個清楚的
模式：**驗證文化非常扎實**（幾乎每支檔案都有具名的 `run_checks()`、雙路徑交叉驗證、對照論文原始
LaTeX 原始碼而非 PDF 或記憶、誠實標記尚未解決的問題），但**協變變分引擎本身（模組 B，長期目標
真正的核心）幾乎沒有進度**——目前的背景/微擾結果全部是繞過協變引擎、直接在 FLRW 背景下用 SymPy
做分量計算得出的。這是合理且誠實的過渡策略，但與長期目標（任意 Lagrangian 自動協變變分）的
方向不同，需要在本規格中明確排出銜接路徑。

另外，目前的背景理論工作對照的是 arXiv:1603.05806（較早、較簡單版本的廣義 Proca：G2(X) 不含
F、Y；L5 沒有 g5(X) 項），而不是你的近期驗收目標 arXiv:1703.09573（G2(X,F,Y)；L5 多一項）。
兩篇論文作者重疊、理論有承接關係，但作用量不是同一個。本規格把「升級背景理論到 1703.09573 的
真正版本」列為階段一必須完成的工作，而非可以延用既有結果直接宣稱完成。

### 目標

**長期目標（不可量化，但要有可量化的階段性代理指標）**：建立一套能對任意由度規、黎曼張量及其
高階不變量、協變導數、以及額外物質場構成的重力作用量，自動做變分、導出修改後場方程式的程式
系統，使用者只需要寫下 Lagrangian 密度即可，不需要人工逐項手推。

**近期驗收目標（可量化，見下方「驗收條件」章節，對應你原規格的五個階段）**：
用本系統重現 arXiv:1703.09573v2 第 II、III 節的全部指定方程式（背景 (2.11)-(2.13)；張量微擾
(3.3)-(3.5)；向量微擾 (3.7)-(3.8)；純量微擾 (3.9)-(3.27)），每一條都必須同時滿足：
1. 對照到論文哪一條方程式編號（不是自己重新定義一個等價但無法逐項比對的形式）；
2. 有兩條獨立推導路徑互相驗證（見「H2 雙獨立路徑」定義）；
3. 沒有任何步驟是「抄公式」——凡尚未獨立驗證的文獻公式，在系統與本規格中都必須標記為
   `provenance: literature_unverified`，不得計入「已完成」。

### 範圍

**本次規格涵蓋（全部 A–J 模組，五個驗收階段全部在範圍內）**：
- 模組 A：全案唯一的慣例/理論定義層。
- 模組 B：協變變分引擎（δΓ 原語、自動 IBP、曲率不變量變分規則庫、G_i(X,F,Y) 鏈式法則自動化）。
- 模組 C：不變量建構庫（X, F, Y, F̃, L^{μναβ}, G_{μν} 從 g, A, ∇ 協變建構）。
- 模組 D：cadabra 安全層（已有雛形，擴大套用範圍到全部 cadabra 腳本）。
- 模組 E：背景層（升級到 1703.09573 真正版本，補齊 L2 的 F, Y 依賴與 L5 的 g5 項）。
- 模組 F：微擾層（統一的 ansatz 建構器、SVT 投影算子庫、符號化二階展開機制、通用約束求解、
  通用 N×N kinetic matrix 萃取）。
- 模組 G：真正的 Schutz-Sorkin 物質部門（純量部門要用真正的 J^μ 推導取代目前的 k-essence 對偶）。
- 模組 H：驗證架構（把現有零散的 `run_checks()` 收斂成可一鍵執行、有回歸保護的測試網）。
- 模組 I：可追溯性（docstring 慣例已經不錯，補上 LaTeX 輸出能力）。
- 模組 J：效能/擴展性（在 B 模組真正做出協變引擎後，逐步把高風險的分量計算換成抽象指標）。

**不做什麼**：
- 不做圖形介面、不做 Web 服務、不做多人協作、不做雲端部署。這是本機/WSL 執行的研究程式庫。
- 不追求對任意論文、任意理論的通用重現能力——本規格的近期驗收範圍限定在 arXiv:1703.09573
  第 II、III 節；階段五的「全新理論」測試只需要一個示範案例（如 Gauss-Bonnet 或 Riemann²），
  不追求覆蓋所有可能的高階曲率不變量。
- 不做效能極致優化；J 模組的要求是「不要因為選錯方法而組合爆炸」，不是「要多快」。
- 不引入資料庫、不引入網路服務依賴；所有「持久化」都是 git 版控的原始碼與登記檔（見下方
  「核心資料模型」）。

### Persona

本專案只有兩類使用者，兩者都需要被規格明確服務：

| Persona | 角色 | 主要需求 | 對規格的具體要求 |
|---|---|---|---|
| 你（物理教授，專案負責人） | 唯一的人類決策者與最終審核者 | 快速判斷某個推導是否可信、卡在哪、要不要接受某個 workaround | 每個模組的輸出都要能讓你在不重新推導的情況下，看出「這是獨立驗證過的」還是「這是抄文獻的」；卡住時要看到「試過什麼、為什麼判斷走不通」，不是只看到失敗訊息 |
| 程式碼協作 agent（Claude Code / Codex，跨對話 session） | 實際寫程式、跑腳本、回報結果的執行者 | 需要在沒有先前對話記憶的情況下，快速掌握慣例（模組 A）、知道哪些是已驗證的地基、哪些是禁區（例如不能呼叫 `rename_dummies()`） | 每個模組要有清楚、自足的 docstring 與慣例文件；D 模組的三條規則必須是「程式介面」而不是「註解裡的提醒」（你的原規格原文），本規格延續這個要求並擴大到 A、B、C、G 模組的慣例 |

### 系統說明

系統分成五層，由下往上是嚴格依賴關係（上層只能透過下層公開介面存取，不得繞過）：

```
┌─────────────────────────────────────────────────────────────┐
│ 第 5 層 驗證與追溯（模組 H、I）——橫向貫穿全部層級，不屬於依賴鏈  │
├─────────────────────────────────────────────────────────────┤
│ 第 4 層 物質部門（模組 G）：Schutz-Sorkin 流體、Proca 向量場    │
├─────────────────────────────────────────────────────────────┤
│ 第 3 層 應用層：背景 mini-superspace（E）／微擾層（F）           │
│         兩者都是「協變場方程 + 特定 ansatz」的代入與化簡         │
├─────────────────────────────────────────────────────────────┤
│ 第 2 層 協變變分引擎（模組 B）+ 不變量建構庫（模組 C）           │
│         B 消費 C 建出的不變量，對它們做 δ(√-g L)/δg, δ(場)     │
├─────────────────────────────────────────────────────────────┤
│ 第 1 層 慣例/理論定義層（模組 A）+ cadabra 安全層（模組 D）      │
│         全案唯一符號來源；D 是操作 cadabra2 的強制介面          │
└─────────────────────────────────────────────────────────────┘
```

第 1 層是地基：模組 A 定義「這個專案裡 X 是什麼、正負號怎麼訂」，模組 D 定義「怎麼安全地叫
cadabra2 而不撞上它的三個已知 bug」。第 2 層是長期目標的核心：模組 C 提供「用 g, A, ∇ 協變
建構出 X, F, Y, F̃, L^{μναβ}, G_{μν}」的建構器，模組 B 提供「對含有這些不變量的任意純量
Lagrangian 做變分」的引擎，兩者合起來就是「餵 Lagrangian 進去、吐場方程出來」的長期目標系統。
第 3 層是把第 2 層的協變結果代入具體 ansatz（FLRW 背景、SVT 微擾分解）得到可與論文比對的具體
方程式——這一層目前的實作（sympy_layer 全部 12 支檔案）繞過了第 2 層直接用分量計算，是本規格
要收斂、但不要求立刻整層重寫的過渡狀態（見「開發應注意重點」）。第 4 層是把 Schutz-Sorkin
流體、Proca 向量場接進第 3 層的 ansatz。第 5 層橫向貫穿：任何一層的輸出都要能被 H 模組拿去跟
論文比對、被 I 模組標上來源。

### 核心流程設計

**流程一：新增一個 Lagrangian 項並取得場方程（長期目標的核心使用情境）**

1. 使用者（或協作 agent）在模組 A 的理論登記檔（`conventions/gp_conventions.py`）裡，用模組 C
   的建構器語彙寫下新的 L_i（例如 `L6 = G6(X) * dual_riemann(mu,nu,al,be) * cov_d(A,nu,mu) * cov_d(A,be,al)`）。
2. 呼叫模組 B 的 `vary_metric(L_i)` 與（若含物質場）`vary_field(L_i, A)`：
   - 內部依序執行：分配全新啞指標配額（D3）→ 對每個不變量因子套用 B4 規則庫的已知變分規則
     （若是複合不變量，先用 C 模組的建構規則展開到 g, A, ∇ 層級）→ 若含巢狀 ∇，套用 B2 的
     δΓ 原語 → B3 自動部分積分並列出被丟棄的邊界項 → 用 D1/D2 規則安全地把各項相加、canonicalise。
   - 每一步都印出中間式，讓使用者可以在任何一步中斷檢查（呼應「我對工作方式的要求」第 2 點）。
3. 輸出：{場方程式（cadabra Ex 物件）、被丟棄的邊界項列表、每個不變量因子用到的變分規則來源（B4
   規則庫的哪一條）、L_i 的 provenance 標記}。
4. 若要代入具體 ansatz 驗證（第 3 層），呼叫 E 或 F 層的 `substitute_ansatz(field_eq, ansatz)`，
   同時獨立跑一次「直接對 mini-superspace/二階作用量做 Euler-Lagrange」（H2 的第二條路徑），
   兩者必須符號上完全一致，否則整個流程回報 `MISMATCH`，不得静默採信任一邊。

**流程二：重現論文一條方程式（近期驗收目標的核心使用情境）**

1. 在 `verification/paper_equations/arxiv_1703_09573.py` 裡先把目標方程式逐字轉錄成符號表達式
   （H1，論文即標準答案），標記來源頁碼/式號、轉錄依據（LaTeX 原始碼 or 官方勘誤，不接受 OCR
   或記憶）。
2. 執行流程一得到系統自己推導出的表達式。
3. 呼叫 H 模組的 `assert_symbolically_equal(derived, paper_ground_truth, tolerance=...)`；數值
   抽樣交叉檢查（H4）作為符號比對的輔助證據，不是唯一證據。
4. 若通過，登記進 `verification/registry.json`（模組 I 的來源可追溯清單），標記
   `{eq_ref, derivation_route_A, derivation_route_B, verified_date, status: reproduced}`。
5. 若不通過：先假設是自己錯，依你規格「四」的要求徹底排查；仍不符時，把兩邊推導都攤開輸出
   成報告（不是只印 diff），交給你判斷。

### 開發應注意重點以及應避開誤區

以下每一條都直接對應既有程式碼裡已經真實發生過、且被誠實記錄下來的教訓，收斂進本規格作為
強制規則，避免未來重蹈覆轍：

1. **絕不在各檔案裡各自定義符號慣例。** 現有 `proca_minisuperspace.py`、`adm_scalar.py`、
   `schutz_sorkin_vector.py`、`fluid_sector.py`、`proca_tensor.py`、`vector_field_equation.py`
   至少 6 支檔案各自重新宣告 `t`, `k`, `Mpl`, `X` 等符號——這正是你曾經歷的「G5,X 正負號不一致」
   的根源模式。模組 A 落地後，這些宣告一律改成 `from conventions.gp_conventions import t, k, Mpl, X_of`。
2. **cadabra 安全層三條規則必須套用到每一支新 cadabra 腳本，沒有例外**（延續現有
   `cadabra_utils.py` 的規則，見模組 D）：不讓兩個帶 `\delta{...}` 的項同時出現在同一個
   `Ex()`/`substitute()` 字串裡；只用 `canonicalise()`，禁止單獨呼叫 `rename_dummies()`；新
   dummy pair 一律從 `fresh_indices()` 取名。`check_no_rename_dummies.py` 這類 lint script
   要在每次新增 cadabra 檔案後執行，且應該收進 Stage 0 的自動檢查（模組 H）。
3. **不接受用「比較好算但不是原問題的模型」代替真正的物理對象。** 目前
   `scalar_plus_fluid_sector.py` 用 k-essence 對偶取代真正的 Schutz-Sorkin 純量流體，這是本規格
   明確判定違反你規格 G 段要求、必須重做的項目，不能因為「數值上接近」而繼續沿用。
   `schutz_sorkin_vector.py` 的向量部門做法（真正從 J^μ 出發，逐項比對論文式子）才是正確範本。
4. **文獻抄錄的公式必須誠實標記，不得混入「已驗證」。** `proca_tensor.py` 對 G5 相關項的自白
   （"G5-dependent terms' FORM ... is still quoted from the literature rather than independently
   re-derived"）是好範本，應該推廣成模組 I 的強制 docstring 欄位 `provenance`，而不是個案自覺。
5. **L4/L5 的 g^{μν} 變分需要完整 δΓ 展開，不能用「已縮併形式」的技巧繞過。**
   `proca_L3.py` 自己的檔尾註記已經指出：f(R) 與 L3 的 A_μ 變分能繞過 δΓ，是因為被變分對象本身
   不被 ∇ 微分；L4（兩層巢狀）、L5（三層巢狀）沒有這個捷徑，必須先把 B2 的 δΓ 原語做穩，才能
   繼續往上蓋。這是判斷「現在該做 B2 還是該先湊近期驗收」時的關鍵技術事實。
6. **IBP 不是「試一次」就一定收斂。** `proca_minisuperspace.py` 的 `ibp_remove_addot()` 文件
   已經記錄過一次真實的死路（G5(X)G_{μν}∇^μA^ν 項用係數提取法怎麼迭代都不收斂，換用 Bianchi
   恆等式才解決）。B3 的自動 IBP 機制設計時，必須把「提取係數法不保證收斂，需要有替代識別式
   （如 Bianchi 恆等式）的介面」當成一等公民案例，不能只支援最簡單的單一 IBP pass。
7. **線性化不能只做到一階就平方。** `adm_scalar.py` 的三個 bug 修復紀錄（K 的線性化、X 的
   g^{00} 展開漏掉二階項、abar-chibar 交叉項升冪錯誤）全部源自「先線性化再平方/相乘」這個
   捷徑。F3（符號化二階展開機制）必須是「先展開到二階、再平方/相乘」，且要有 H4 數值抽樣交叉
   檢查作為預設步驟，不是事後才想到要做。
8. **不要用 `sp.series()` 對多變數依序展開。** `schutz_sorkin_vector.py` 的文件記錄過
   `sp.series()` 逐變數截斷會漏掉聯合階數超過但個別變數階數未超過的交叉項（如 V_i²W_i²）。
   F3 一律採用「顯式 eps 標記 + 對 eps 取 Taylor 係數」的模式，這是目前程式碼裡唯一被驗證過
   安全的展開方法，應該收進模組 F 的標準介面而不是各檔案各自重新發明。
9. **不要假設環境已裝好。** cadabra2 只在 WSL Ubuntu 可用（Windows 端的 `py -3` 沒有它），
   sympy 只在 Windows 端的 `py -3` 可用（本次調查時系統預設的 `python3` 兩者都沒有）。任何
   Codex/Claude Code 階段計畫的「測試要求」都必須寫明要在哪個直譯器下跑，不能假設單一 python
   指令兩層都能執行。

### 執行環境

| 項目 | 版本/位置 | 驗證方式 | 備註 |
|---|---|---|---|
| cadabra2 | 2.5.14，WSL Ubuntu 系統套件（`/usr/lib/python3.12/dist-packages/cadabra2...so`） | `wsl.exe -- python3 -c "import cadabra2"` 實測成功；`dpkg -l` 確認版本 | 目前是 WSL 開機後才啟動（Stopped 狀態，隨呼叫自動啟動）；不建議切換到其他安裝方式，因為 D 模組三條 workaround 是針對這個版本重新驗證過的（2026-09-12） |
| Python + SymPy | Windows `py -3`，sympy 1.14.0 | `py -3 -c "import sympy"` 實測成功 | 系統預設 `python3`／`python` 兩者都**沒有** sympy，也沒有 cadabra2——任何腳本呼叫範例都必須明確寫 `py -3 script.py`，不能只寫 `python3 script.py` |
| 版本升級政策 | 見上方「cadabra2 版本調查結果」 | — | 只有在上游 changelog 明確聲稱修復過 D 模組列出的三個具體行為時才升級；否則固定在 2.5.14，避免未經驗證的版本漂移 |

### 任務模型與資訊優先級（推導任務版，非 UI 版）

本專案沒有互動式工作台介面，因此不適用 output-template 原版的「首屏／視覺群組／主 CTA」框架。
以下改用「推導任務」的等價結構，服務相同的目的：讓任何一次協作 session（人類或 agent）能立刻
判斷現在該優先處理什麼、什麼可以先擱置。

#### 任務模型表

| 層級 | 內容 | 為何屬於這一層 | 是否為目前 session 的優先項 |
|---|---|---|---|
| 唯一主目標 | 完成階段一：把 E 模組背景理論升級到 1703.09573 真正版本（G2(X,F,Y)、補 g5 項），逐項比對 (2.11)-(2.13) | 是全案地基，其他階段（張量/向量/純量微擾）都要在正確的背景解上展開，背景理論選錯論文版本會讓後面全部重做 | 是（見下方建議執行順序第 2 點） |
| 次目標 | 模組 A 落地（慣例單一來源）、模組 D 擴大套用範圍、cadabra B2（δΓ 原語）在 L3 的 g^{μν} 變分上驗證 | 降低未來重工風險，且 B2 是長期目標的核心，越早驗證越早知道可行性 | 是，與主目標並行 |
| 低頻目標 | G 模組純量流體重做（真 Schutz-Sorkin）、F 模組三部門收斂成單一 ansatz 建構器 | 需要等主目標的正確背景理論定案後才有穩定的地基可以站 | 否，排在主目標之後 |
| 罕見目標 | 階段五：套用到全新理論（Gauss-Bonnet/Riemann²）、J 模組把分量計算換成抽象指標 | 是長期目標的最終驗收，但前置條件（B 模組完整）尚未成立，現在做風險是重工 | 否，明確排在最後 |

#### 資訊分類表

| 資訊項目 | 分類 | 使用頻率 | 是否必須立刻可見 | 不顯示的風險 |
|---|---|---|---|---|
| `run_checks()` PASS/FAIL 清單 | status-feedback | 每次跑腳本都要看 | 是 | 看不到會誤以為全部通過，重蹈"6/7 通過但沒人注意那 1 項"的風險 |
| 已知未解決問題（如 A0 方程式殘留 addot 項） | exception-handling | 低頻，但一旦存在必須顯眼 | 是（不可被其他 PASS 訊息淹沒） | 被忽略會讓下一位協作者誤用一個已知有問題的方程式當地基 |
| provenance 標記（derived / literature_unverified） | decision-supporting | 每次要判斷"能不能算完成"時查看 | 是 | 混入「已驗證」清單會讓驗收階段的判斷失真 |
| D 模組三條 cadabra workaround 規則 | reference | 每次寫新 cadabra 腳本時查一次 | 否，但要容易找到（建議放在 `cadabra_utils.py` 頂部 docstring，現況已經如此） | 忘記查會直接撞上已知 bug，重複踩雷 |
| cadabra2 版本調查歷史紀錄 | audit-history | 極低頻，只有考慮升級版本時查 | 否 | 不影響日常工作，但升級決策時若遺失會重做整個調查 |

### 狀態模型與揭露策略（推導任務狀態機，非 UI 版）

每一條「論文方程式重現任務」都走以下狀態機（取代 UI 版的 state model）：

#### 狀態矩陣

| State | 進入條件 | 這個狀態代表什麼 | 必顯資訊 | 允許暫時隱藏 | 離開條件 |
|---|---|---|---|---|---|
| `not_started` | 尚未動手 | 論文式子已登記進 `paper_equations/`，但還沒推導 | 對應的論文式號、頁碼 | 推導細節 | 開始推導 |
| `derived_single_route` | 完成一條推導路徑 | 有結果，但只有一條推導路徑，尚不可信 | provenance=derived，但需標註「僅一條路徑」 | 與論文的逐項比對結果 | 補上第二條獨立路徑 |
| `self_consistent` | 兩條獨立路徑（H2）互相吻合 | 內部自洽，但還沒跟論文比對過 | 兩條路徑各自的推導摘要 | — | 與論文 ground truth 比對 |
| `paper_matched` | 通過 H1（對照論文逐項比對）與 H4（數值抽樣）雙重檢查 | 可以宣稱「重現成功」（符合你規格三.3 的要求） | 對照的論文式號、兩條獨立路徑名稱、比對方式 | — | 進入回歸測試集 |
| `regression_locked` | 已登記進 `verification/registry.json`，且被 H 模組的一鍵 runner 涵蓋 | 任何未來修改若破壞這條結果，會被自動攔截 | 上次驗證日期、涵蓋它的測試檔案路徑 | — | （終態，除非上游 cadabra2 版本改變才重新開） |
| `known_issue` | 檢查跑過但目前未通過（如 A0 方程式殘留 addot） | 明確、誠實地卡住，不是靜默失敗 | 具體卡在哪一步、試過什麼、目前假設的可能原因 | — | 問題解決後轉回 `derived_single_route` 重新走一次 |

#### 揭露策略（取代 UI 版的容器/收合規則）

- `run_all_checks.py`（模組 H 新增的一鍵 runner）的預設輸出，只顯示每條任務目前的 state 與一行
  摘要（對應「必顯資訊」欄）；推導細節、兩條路徑各自的完整符號式，收在 `--verbose` 旗標或個別
  檔案的 `if __name__ == '__main__'` 輸出裡，避免一次性洪水式輸出把 `known_issue` 淹沒。
- `known_issue` 狀態的任務，在彙總報告中永遠排在最上方，不與 `regression_locked` 的任務混排——
  這是直接回應「6/7 通過但那 1 項不能被忽略」的既有教訓。
- provenance=`literature_unverified` 的項目，在任何彙總報告與 LaTeX 輸出中都要有視覺上不會被
  忽略的標記（例如固定前綴 `[UNVERIFIED]`），不能只在原始碼註解裡出現。

### 專案目錄規劃

```
auto_variation/
├── docs/
│   └── SPEC.md                       # 本規格，全案唯一的範圍/驗收/模組邊界權威來源
├── conventions/                      # 【新增】模組 A：慣例與理論定義層，全案唯一符號來源
│   └── gp_conventions.py             #   signature、Riemann/Ricci 正負號、X/F/Y 定義、
│                                      #   A^mu/A_mu 升降指標慣例、L2..L6 的 G_i(X,F,Y) 登記
├── cadabra/                          # 模組 B（協變半邊）+ 模組 D（安全層）
│   ├── cadabra_utils.py              #   既有：D1-D3 workaround 介面（fresh_indices/sub_copy/…）
│   ├── check_no_rename_dummies.py    #   既有：lint script
│   ├── variation_engine.py           #   【新增】B1-B3：vary_metric()/vary_field()/自動 IBP
│   ├── curvature_invariant_rules.py  #   【新增】B4：δ√-g/δR/δR_munu/δR_riemann/δG_munu 規則庫
│   ├── invariants.py                 #   【新增】模組 C：X/F/Y/F̃/L^{munualbe}/G_munu 建構器
│   ├── eh_variation.py               #   既有，回歸測試用
│   ├── fR_variation.py               #   既有，回歸測試用
│   ├── proca_L2.py / proca_L3.py     #   既有，L3 的 g^{mu nu} 變分待補（見開發重點第 5 條）
│   └── proca_L4.py / proca_L5.py / proca_L6.py   # 【新增】
├── sympy_layer/                      # 模組 E（背景）+ 模組 F（微擾）的分量層驗證與 ansatz 代入
│   ├── tensor_utils.py               #   既有：通用 Christoffel/Riemann/Ricci
│   ├── proca_minisuperspace.py       #   既有，需升級成 1703.09573 版本（G2(X,F,Y)、g5 項）
│   ├── proca_L6_background.py        #   既有：L6 背景貢獻為零的兩個獨立引理證明
│   ├── proca_tensor.py               #   既有：G_T/F_T/c_T^2，G5 項需從文獻抄錄改為獨立重新推導
│   ├── fR_gravity.py                 #   既有：f(R) 雙路徑交叉驗證範本（H2 的優良示範）
│   ├── vector_sector.py              #   既有：只做到 GR+Maxwell，需補 G4X/G5X 的 Z_i 耦合
│   ├── vector_field_equation.py      #   既有：GR 部分的 (3.26)
│   ├── adm_scalar.py                 #   既有：ADM 純量部門地基（三個 bug 修復紀錄）
│   ├── scalar_sector.py              #   既有：GR+Proca+Maxwell 純量部門
│   ├── scalar_plus_fluid_sector.py   #   既有，需整支替換為真正 Schutz-Sorkin（不可沿用 k-essence 對偶）
│   ├── schutz_sorkin_scalar.py       #   【新增】仿照 schutz_sorkin_vector.py 的模式做純量部門
│   ├── schutz_sorkin_vector.py       #   既有：真正 Schutz-Sorkin 範本，比對 (3.20)-(3.23) 過關
│   ├── fluid_sector.py               #   既有：流體背景層
│   ├── pert_engine.py                #   既有：單場零件，需擴充成 F1/F2/F6 的通用版本
│   ├── svt_builder.py                #   【新增】F1：一次產生 (3.1)-(3.2) 全部微擾變數的建構器
│   └── kinetic_matrix.py             #   【新增】F6：通用 N×N kinetic/gradient matrix 萃取
├── verification/                     # 【新增】模組 H：驗證架構
│   ├── paper_equations/
│   │   ├── arxiv_1703_09573.py       #   逐字轉錄的論文 ground truth（H1）
│   │   └── arxiv_1603_05806.py       #   背景理論的舊版對照（保留供比較，不再是驗收依據）
│   ├── run_all_checks.py             #   一鍵 runner，彙總 cadabra/ + sympy_layer/ 全部 run_checks()
│   └── known_limits.py               #   H3：GR/標準 Proca/f(R)/Maxwell/Horndeski 極限測試集中管理
├── traceability/                     # 【新增】模組 I
│   └── latex_export.py               #   把已驗證的 Ex/SymPy 表達式輸出成可直接貼進論文的 LaTeX
├── tests/                            #   既有但空；改成 pytest 包住 verification/run_all_checks.py
│   └── test_registry_regression.py   #   【新增】對 verification/registry.json 裡的 regression_locked 項目跑回歸
├── 1603.05806v2.pdf                  #   既有：背景理論來源論文（供人工核對，不是程式輸入）
├── 1703.09573v2.pdf                  #   既有：近期驗收目標論文
└── .claude/                          #   既有：session 設定，與本規格無直接耦合
```

**命名原則**：cadabra 檔案一律用 `<理論項目>_<動作>.py`（如 `proca_L4.py`）；sympy_layer 檔案
一律用 `<物理部門>_<層級>.py`（如 `adm_scalar.py`, `scalar_sector.py`）；`verification/` 內的
論文轉錄檔一律以 `arxiv_<id去掉點與v版本>.py` 命名，避免像既有程式碼曾發生過的「引用錯論文
式號」問題（`schutz_sorkin_vector.py` 文件記載的 (3.35)-(3.38) 誤植為 (3.20)-(3.23) 事件）。

**測試/資源/設定放置位置**：測試放 `tests/`（新增，目前為空）；沒有執行期設定檔（無環境變數、
無 config.yaml）——執行環境差異完全靠「執行環境」章節裡的直譯器選擇處理，不做成可設定項，
因為兩個直譯器（Windows py -3 / WSL cadabra2）之間不存在「切換」的情境，只有「各自固定用途」。

### 前後端模組（本專案無前後端之分，改為「協變層／應用層」模組邊界）

#### 模組架構圖（SVG）

```svg
<svg viewBox="0 0 1200 760" xmlns="http://www.w3.org/2000/svg">
  <rect x="0" y="0" width="1200" height="760" fill="#ffffff"/>
  <text x="600" y="30" text-anchor="middle" font-size="20" font-weight="bold" fill="#111">廣義重力變分系統 — 模組架構</text>

  <!-- Layer 5: H/I 橫向 -->
  <rect x="30" y="55" width="1140" height="60" rx="8" fill="#fef3c7" stroke="#d97706" stroke-width="1.5"/>
  <text x="600" y="80" text-anchor="middle" font-size="15" font-weight="bold" fill="#92400e">第5層　驗證與追溯（模組 H 驗證架構 / 模組 I 可追溯性）— 橫向貫穿下方全部層級</text>
  <text x="600" y="100" text-anchor="middle" font-size="12" fill="#92400e">verification/run_all_checks.py · paper_equations/ · traceability/latex_export.py</text>

  <!-- Layer 4: G matter -->
  <rect x="30" y="135" width="1140" height="70" rx="8" fill="#dbeafe" stroke="#2563eb" stroke-width="1.5"/>
  <text x="600" y="160" text-anchor="middle" font-size="15" font-weight="bold" fill="#1e3a8a">第4層　物質部門（模組 G）</text>
  <text x="600" y="182" text-anchor="middle" font-size="12" fill="#1e3a8a">schutz_sorkin_vector.py（已驗證）· schutz_sorkin_scalar.py（新增，取代 k-essence 對偶）· fluid_sector.py</text>

  <!-- Layer 3: E/F -->
  <rect x="30" y="235" width="555" height="120" rx="8" fill="#dcfce7" stroke="#16a34a" stroke-width="1.5"/>
  <text x="307" y="260" text-anchor="middle" font-size="15" font-weight="bold" fill="#14532d">第3層a　背景 mini-superspace（模組 E）</text>
  <text x="307" y="282" text-anchor="middle" font-size="12" fill="#14532d">proca_minisuperspace.py（升級至 1703.09573）</text>
  <text x="307" y="300" text-anchor="middle" font-size="12" fill="#14532d">proca_L6_background.py · fR_gravity.py（雙路徑範本）</text>
  <text x="307" y="318" text-anchor="middle" font-size="12" fill="#14532d">目標：與論文 (2.11)-(2.13) 逐項比對</text>

  <rect x="615" y="235" width="555" height="120" rx="8" fill="#dcfce7" stroke="#16a34a" stroke-width="1.5"/>
  <text x="892" y="260" text-anchor="middle" font-size="15" font-weight="bold" fill="#14532d">第3層b　微擾層（模組 F）</text>
  <text x="892" y="282" text-anchor="middle" font-size="12" fill="#14532d">svt_builder.py（新增）· kinetic_matrix.py（新增）</text>
  <text x="892" y="300" text-anchor="middle" font-size="12" fill="#14532d">pert_engine.py · vector_sector.py · adm_scalar.py · scalar_sector.py</text>
  <text x="892" y="318" text-anchor="middle" font-size="12" fill="#14532d">目標：(3.3)-(3.27) 逐項比對</text>

  <!-- Layer 2: B + C -->
  <rect x="30" y="385" width="555" height="110" rx="8" fill="#ede9fe" stroke="#7c3aed" stroke-width="1.5"/>
  <text x="307" y="410" text-anchor="middle" font-size="15" font-weight="bold" fill="#4c1d95">第2層a　協變變分引擎（模組 B）— 長期目標核心</text>
  <text x="307" y="432" text-anchor="middle" font-size="12" fill="#4c1d95">variation_engine.py：vary_metric() / vary_field()</text>
  <text x="307" y="450" text-anchor="middle" font-size="12" fill="#4c1d95">curvature_invariant_rules.py：δΓ 原語、B4 規則庫</text>
  <text x="307" y="468" text-anchor="middle" font-size="12" fill="#4c1d95">現況：僅 EH/f(R)/L2/L3(A變分) 驗證過</text>

  <rect x="615" y="385" width="555" height="110" rx="8" fill="#ede9fe" stroke="#7c3aed" stroke-width="1.5"/>
  <text x="892" y="410" text-anchor="middle" font-size="15" font-weight="bold" fill="#4c1d95">第2層b　不變量建構庫（模組 C）</text>
  <text x="892" y="432" text-anchor="middle" font-size="12" fill="#4c1d95">invariants.py：X, F, Y, F̃, L^munualbe, G_munu</text>
  <text x="892" y="450" text-anchor="middle" font-size="12" fill="#4c1d95">從 g, A, ∇ 協變建構，供模組 B 消費</text>
  <text x="892" y="468" text-anchor="middle" font-size="12" fill="#4c1d95">現況：尚未建立，目前各檔案手刻分量公式</text>

  <!-- Layer 1: A + D -->
  <rect x="30" y="535" width="555" height="100" rx="8" fill="#fee2e2" stroke="#dc2626" stroke-width="1.5"/>
  <text x="307" y="560" text-anchor="middle" font-size="15" font-weight="bold" fill="#7f1d1d">第1層a　慣例/理論定義層（模組 A）</text>
  <text x="307" y="582" text-anchor="middle" font-size="12" fill="#7f1d1d">gp_conventions.py：全案唯一符號來源</text>
  <text x="307" y="600" text-anchor="middle" font-size="12" fill="#7f1d1d">現況：不存在，6+ 檔案各自重複宣告 t/k/Mpl/X</text>

  <rect x="615" y="535" width="555" height="100" rx="8" fill="#fee2e2" stroke="#dc2626" stroke-width="1.5"/>
  <text x="892" y="560" text-anchor="middle" font-size="15" font-weight="bold" fill="#7f1d1d">第1層b　cadabra 安全層（模組 D）</text>
  <text x="892" y="582" text-anchor="middle" font-size="12" fill="#7f1d1d">cadabra_utils.py：D1-D3，已驗證對 cadabra2 2.5.14 有效</text>
  <text x="892" y="600" text-anchor="middle" font-size="12" fill="#7f1d1d">現況：全案完成度最高的模組</text>

  <!-- arrows -->
  <path d="M600 115 L600 135" stroke="#6b7280" stroke-width="2" marker-end="url(#arrow)"/>
  <path d="M600 205 L600 235" stroke="#6b7280" stroke-width="2" marker-end="url(#arrow)"/>
  <path d="M600 355 L600 385" stroke="#6b7280" stroke-width="2" marker-end="url(#arrow)"/>
  <path d="M600 495 L600 535" stroke="#6b7280" stroke-width="2" marker-end="url(#arrow)"/>
  <defs>
    <marker id="arrow" markerWidth="10" markerHeight="10" refX="5" refY="5" orient="auto">
      <path d="M0,0 L10,5 L0,10 z" fill="#6b7280"/>
    </marker>
  </defs>
  <text x="1150" y="740" text-anchor="end" font-size="11" fill="#6b7280">箭頭方向＝依賴方向（上層依賴下層公開介面）</text>
</svg>
```

### 使用流程

見上方「核心流程設計」兩條主流程（新增 Lagrangian 項 / 重現論文方程式）。第三條輔助流程：

**流程三：升級 cadabra2 版本時的重新驗證**
1. 只重跑 `cadabra_utils.py` 文件開頭建議的三個最小重現腳本（Accent 邊界、rename_dummies、
   substitute 命名衝突），不需要重新驗證全部程式碼。
2. 若三個 bug 有任一個行為改變（修復或惡化），更新 `cadabra_utils.py` 的 docstring 與「執行
   環境」章節的版本調查紀錄，並評估是否需要調整 D 模組的 workaround 規則。
3. 若三個 bug 行為不變，只更新「上次重新驗證日期」，不需要改動任何 workaround 邏輯。

### 功能清單（含物件生命週期，本領域的 CRUD 等價物）

本專案沒有使用者建立的「物件」在傳統 CRUD 意義上，但有三類需要明確生命週期管理的登記物件：

| 物件 | 新增（Create） | 更新（Update） | 棄用/移除（Delete） | 狀態欄位 |
|---|---|---|---|---|
| B4 曲率不變量變分規則（`curvature_invariant_rules.py` 裡的一條規則） | 新增一個高階不變量時，只需要登記它的變分規則（呼應你規格 B4 原文「不是重寫整條推導鏈」） | 若發現既有規則有誤（如曾發生過的 G5,X 正負號事件），更新規則並觸發全部依賴它的 regression 重跑 | 規則一旦被論文比對驗證過，不應移除，只能被更精確的推導取代（保留舊版於 git 歷史） | `provenance`（derived / literature_unverified）、`verified_against`（論文式號或 None） |
| D 模組 cadabra bug workaround 條目 | 發現新的 cadabra2 行為異常時新增 | cadabra2 版本升級後，若行為改變需更新條目與重新驗證日期 | 只有在官方 changelog 明確聲稱修復後才移除，且移除前要重新驗證 | `bug_id`、`last_reverified_version`、`last_reverified_date` |
| `verification/registry.json` 的方程式重現條目 | 一條方程式通過流程二第 4 步後新增 | 上游模組（A/B/C/E/F/G）有變動導致重新推導時更新 | 不刪除，只能狀態轉移（見狀態矩陣，`known_issue` 不等於刪除） | 見上方「狀態模型」章節的完整狀態機 |

### 非功能需求

| 面向 | 要求 | 對應既有教訓 |
|---|---|---|
| 效能（模組 J） | cadabra 腳本單次執行時間應在數分鐘內完成（合理上限，非硬性 SLA）；若巢狀 ∇ 展開導致組合爆炸，優先檢討是否該用抽象指標路線而非加大分量計算規模 | 現有策略「用 sympy 分量計算頂著」在 FLRW 背景下可行，但不可無限期擴大到非對角/含空間依賴的微擾度規 |
| 可靠性 | 任何 `run_checks()` 的 PASS 宣告必須是可重跑、確定性的（不依賴隨機種子除非明確做數值抽樣，且數值抽樣要多組種子交叉驗證） | `fR_gravity.py` 的 `numeric_compare()` 已經是正確範本（多組隨機時間點與係數比對，不依賴單一巧合值） |
| 可維護性 | 新增一個 L_i 項，理論上只需要動 `conventions/gp_conventions.py`（登記）與呼叫模組 B 的通用介面，不應該需要修改 B/C 模組本身的程式碼 | 目前每個 L_i 都要重新手刻整條推導鏈（proca_L2.py, proca_L3.py 各自獨立），違反這條原則，是階段五驗收的核心考驗 |
| 可稽核性（模組 I） | 每個已驗證公式都能一鍵輸出 LaTeX，附上來源與 provenance | 目前完全沒有 LaTeX 輸出能力，是本規格新增的落差 |
| 正確性優先於速度 | 任何「化簡」在採用前必須先被獨立驗證過，不能因為讓後續運算變快而跳過驗證步驟 | `proca_minisuperspace.py` 的「P1=MM_trace、Q1=Q2」化簡是先驗證再採用的正確範例，應作為 J 模組所有化簡的標準流程 |

### 核心資料模型

本專案沒有執行期資料庫；所有「資料」都是 git 版控的 Python 原始碼與人類可讀的登記檔。以下
定義三個核心「實體」，作為 `verification/registry.json`（新增）與各模組 docstring 慣例的欄位規格：

**實體一：`DerivedFormula`（一條已推導或待推導的公式）**
| 欄位 | 型別 | 說明 |
|---|---|---|
| `id` | string | 唯一鍵，慣例 `<paper_id>_eq<eq_number>`，如 `1703_09573_eq3_5` |
| `paper_ref` | string | 對應論文與式號，如 `"arXiv:1703.09573 Eq (3.5)"` |
| `ansatz` | string | 假設的背景/規範，如 `"flat FLRW, A^mu=(A0(t),0,0,0), flat gauge"` |
| `provenance` | enum | `derived` / `literature_unverified`；不得有第三種模糊狀態 |
| `derivation_route_a` | string | 第一條獨立推導路徑的描述（模組/函式路徑） |
| `derivation_route_b` | string \| null | 第二條獨立推導路徑；`self_consistent` 以上狀態必須非 null |
| `state` | enum | 見「狀態模型」章節的六個狀態 |
| `last_verified_date` | date \| null | — |
| `known_issue_note` | string \| null | `known_issue` 狀態時必填：卡在哪、試過什麼 |

**實體二：`ConventionEntry`（模組 A 登記的單一符號慣例）**
| 欄位 | 型別 | 說明 |
|---|---|---|
| `symbol` | string | 如 `X`, `A_mu`, `signature` |
| `definition` | string | 如 `X = -1/2 A_mu A^mu` |
| `sign_note` | string | 明確寫出這個慣例與常見替代慣例的差異，避免未來對照文獻時混淆 |

**實體三：`BugWorkaround`（模組 D 登記的單一 cadabra2 已知行為異常）**
| 欄位 | 型別 | 說明 |
|---|---|---|
| `bug_id` | string | 如 `accent_boundary_free_index_checker` |
| `symptom` | string | 觀察到的錯誤訊息或行為 |
| `workaround_rule` | string | 強制介面規則（如 D1/D2/D3） |
| `cadabra_version_tested` | list[string] | 如 `["2.4.5.4", "2.5.14"]` |
| `last_reverified_date` | date | — |
| `upstream_issue_url` | string \| null | 若已回報上游，填 issue 連結（目前三個 bug 皆為 null，建議依本規格第 0 節開 issue） |

### State 管理與持久化（本領域等價：一切狀態即原始碼 + 登記檔）

本專案沒有執行期 session、沒有使用者草稿需要保存。最接近「可續編進度」的持久化物件是
`verification/registry.json`（新增）：它記錄每條 `DerivedFormula` 目前的 `state`，任何協作
agent 在新 session 開始時，第一步就是讀這個檔案，而不是重新從對話記錄猜測進度。git commit
本身就是持久化機制——每完成一個 state 轉移，對應一次 commit（延續現有 11 筆 commit 的模式，
每筆訊息已經包含類似「6/7 checks pass, 1 flagged open」的狀態摘要，這個慣例應該延續並正規化
成 `registry.json` 的結構化版本，而不是只留在 commit message 的自由文字裡）。

### 錯誤處理／回退策略／可觀測性

**錯誤處理**（符號計算領域的錯誤，取代傳統的例外處理/HTTP 錯誤碼）：

| 錯誤情境 | 處理方式 |
|---|---|
| cadabra2 撞到 D 模組列出的已知 bug | 不是「例外」，是「忘記套用 D1-D3 規則」——`check_no_rename_dummies.py` 這類 lint 在此情境下應該視為錯誤處理的第一道防線，而非事後除錯 |
| H2 雙路徑推導結果不一致 | 不得靜默採信任一邊或取平均；回報 `MISMATCH`，兩條路徑的完整推導都要輸出，交由人工判斷（你規格「四」的要求） |
| IBP 不收斂（如 `ibp_remove_addot` 迭代超過 `max_passes`） | 明確拋出 `RuntimeError`（既有程式碼已經這樣做），不得回傳部分結果假裝成功 |
| 已知未解決問題（如 A0 方程式殘留 addot 項） | 明確標記 `known_issue`，在彙總報告中永遠可見，不得因為其他檢查通過而被稀釋 |

**可觀測性**：`run_all_checks.py`（新增）輸出兩種格式——人類閱讀的終端機文字報告（延續現有
PASS/FAIL 慣例）與機器可讀的 JSON 摘要（供協作 agent 判斷目前整體進度，不需要 regex 解析文字
輸出）。JSON 摘要至少包含：每個 `DerivedFormula` 的 `id`、`state`、上次驗證日期；彙總的
`known_issue` 清單（永遠排最前面）；`literature_unverified` 清單。

### 長時間運算與重試策略（取代 UI 版「通知與背景執行」）

本專案沒有非同步工作佇列，但巢狀 ∇ 展開（尤其 L4/L5/L6 的 δΓ 全展開）有真實的組合爆炸風險
（模組 J 的核心關切）。規則：

- 任何 cadabra 腳本若預期執行時間可能超過幾分鐘（例如 L5 的三層巢狀 δΓ），必須在腳本內把
  推導拆成可獨立重跑的中間步驟（延續現有 `eh_variation.py`/`fR_variation.py` 用編號
  `print('1)', v)` 標記中間結果的慣例），讓超時後可以從中間步驟續跑，而不必整個重來。
- 不設自動重試機制——符號計算若失敗，重跑不會產生不同結果，需要的是修正推導或套用 D 模組
  workaround，自動重試只會浪費時間並可能掩蓋真正的問題。
- 若一個推導在合理時間內（協作 agent 應自行設定不超過 10 分鐘的合理上限）沒有得出結果，
  應中止並回報「卡在哪一步、目前的中間結果是什麼」，而不是無限期等待。

### 建議補充的功能

- `verification/run_all_checks.py` 的 JSON 輸出格式，未來可以直接餵給 `docs/SPEC.md` 自動
  產生一份「目前驗收進度」摘要表，減少人工手動更新規格文件的負擔（本次規格本身尚未做這個
  自動化，列為未來優先度中等的延伸功能）。
- `traceability/latex_export.py` 除了輸出單條公式，未來可以考慮直接輸出成可貼進論文附錄的
  完整推導章節（含每一步驟的中間式），優先度中等，待模組 I 基本功能穩定後再評估。

### 驗收條件

延續你原規格「四、驗收標準」的五個階段，補上本規格新增的可量化細節：

| 階段 | 驗收內容 | 可量化判準 |
|---|---|---|
| 階段一 | **【2026-09-18 完成】** 背景方程 (2.11)-(2.13)（arXiv:1703.09573 be1-be3）由 L2-L5 推導後代入 FLRW 得到並逐項比對論文。（g5、G6、G2 的 F,Y 依賴，依論文原文 astrophv2.tex 第 323-326 行明確聲明對背景無貢獻，見「風險與未決事項」item 2 更新——不需要、也不應該加進背景推導） | `verification/registry.json` 中 `1703_09573_be1/be2/be3` 三條 `DerivedFormula` 皆為 `paper_matched`，H1（逐字比對 `verification/paper_equations/arxiv_1703_09573.py`）與 H2（`sympy_layer/proca_background_direct_route.py` 獨立協變路徑，符號完全相等）雙重驗證通過。過程中發現並修復三個既有 bug：`euler_lagrange_all()` 的 EL_N 漏了 d/dt(∂L/∂Ṅ) 項、EL_a 漏了 Ostrogradsky 二階修正項、L5 的 Bianchi 恆等式化簡捷徑本身不成立（已改用直接 Einstein 張量計算）——原本的「A0 方程式殘留 addot」known issue 因此解決 |
| 階段二 | 【目前狀態：未開始，仍是 GR+Maxwell 子集】(3.3)/(3.4)/(3.5) 張量部門、(3.7)/(3.8) 向量部門，由二階作用量推導得出 | 對應 `DerivedFormula` 達到 `paper_matched`；向量部門必須包含 G4X/G5X 的 Z_i 耦合項（現況 `vector_sector.py` 只做到 GR+Maxwell 子集，不足以宣稱完成）。`sympy_layer/kinetic_matrix.py`（F6）已可用於這個部門一旦二階作用量做出來之後的矩陣萃取，但作用量本身尚未推導 |
| 階段三 | 【目前狀態：未開始】(3.16)-(3.23) 的 w1-w7 以一般 G_i 符號建出，據此導出 (3.10)-(3.15)，並能在標準 Proca 極限下化約 | 對應 `DerivedFormula` 達到 `paper_matched`；標準 Proca 極限測試收進 `verification/known_limits.py`（H3，已建立，該項目標記 `covered=False`，誠實記錄尚未做）作為常設回歸測試 |
| 階段四 | 【2026-09-18 純量流體子部分完成，kinetic matrix 部分仍待 Stage 5-6】(3.24)-(3.27) 的 Q_S、c_S²、q_S、μ_S 由二場 kinetic matrix 自動萃取而得 | 真正 Schutz-Sorkin 純量流體（`sympy_layer/schutz_sorkin_scalar.py`）已從 J^μ 定義從頭推導，對照 arXiv:1605.05066（1703.09573 自己引用的原始出處）的 eq (SMS) 與 eq (deltaj) 逐項通過，取代 `scalar_plus_fluid_sector.py` 的 k-essence 對偶（該檔案已加註淘汰說明但保留供參考）。**誠實記錄一個已知、已理解但尚未解決的殘差**：純量流體矩陣作用量裡有一項純背景（α²、χ²，不耦合任何物質微擾）的殘差，論文引用的 (SMS) 式沒有這項——結構上判斷它屬於重力+物質「合併」二階作用量的 w4 類係數（該合併作用量裡確實有一個「+w4α²」項），需要 GP 重力部門自己的 α,χ 二階作用量才能驗證這個猜測，屬於 Stage 5（F 模組）範圍，見 `verification/registry.json` 的 `1605_05066_SMS_background_residual`（標記為 `known_issue`，未計入已完成）。【2026-09-18 補充】F6 通用 N×N kinetic matrix 萃取器（`sympy_layer/kinetic_matrix.py`）已建立並通過 H2 驗證（對照既有 `scalar_plus_fluid_sector.py` 手刻的 2×2 K_VV/K_SS/K_VS，含非對角交叉項，完全符號相符）；`sympy_layer/svt_builder.py`（F1，統一微擾變數宣告層，採用論文 (3.1)-(3.2) 的確切命名，附新舊變數名對照表）也已建立並通過 H2 驗證。但尚未把 kinetic_matrix.py 套用到「真正 Schutz-Sorkin 純量流體 + GP 重力部門」耦合後的完整二場系統（那需要先把 Stage 4 標記為 known_issue 的重力+物質合併二階作用量做出來），也**還沒**把 vector_sector.py/adm_scalar.py/scalar_sector.py 三支既有檔案改成呼叫 svt_builder.py（這是規格裡明確提醒的高風險大範圍重構，刻意留到之後逐檔案個別驗證，不在這次一次性批次完成） |
| 階段五 | 【目前狀態：未開始，前置條件部分具備】同一套引擎套用在全新高階曲率作用量（如 Gauss-Bonnet 或 Riemann²）上，無人工介入產生修改後場方程 | 選定一個示範案例，從 `conventions/gp_conventions.py` 登記到 `variation_engine.py` 輸出場方程全程不修改 B/C 模組本身程式碼；這是長期目標達成的直接證據。現況：`cadabra/variation_engine.py` 的 δΓ 原語已在 L3（單層巢狀）驗證可行（見階段一下方 Stage 3 記錄），但尚未驗證能否撐住 L4/L5/Gauss-Bonnet 這類更深層巢狀或更高階曲率不變量，`curvature_invariant_rules.py`（B4 規則庫）也還沒建立 |

### 測試案例（節錄，完整清單隨 Stage 2-8 逐步落地於 `tests/`）

- `test_A0_algebraic`: 對任意 G2-G6 選擇，背景 A0 場方程式必須不含 addot 與 Addot（回歸現有
  `proca_minisuperspace.py` 的 check 1，升級後的版本必須保留這條）。
- `test_L6_background_zero`: 對任意 FLRW + 純時間向量 ansatz，L6 對背景的貢獻恆為零（回歸
  `proca_L6_background.py` 的兩個引理）。
- `test_dual_route_consistency`: 對每一條 `state >= self_consistent` 的 `DerivedFormula`，
  兩條獨立推導路徑的符號差必須化簡為零。
- `test_no_rename_dummies_violation`: 對 `cadabra/` 全部檔案跑 `check_no_rename_dummies.py`，
  必須零違規（Stage 0 起就要能跑，並收進 CI 等價的一鍵檢查）。
- `test_schutz_sorkin_scalar_no_kessence_shortcut`: 靜態檢查 `schutz_sorkin_scalar.py` 的
  Lagrangian 是否真的從 `J^mu` 建構（檢查是否出現 `J` 相關符號），防止未來又退回 k-essence 對偶。
- `test_known_limits`: GR、標準 Proca、f(R)、Maxwell 四個極限（Horndeski 對應極限待模組 B 完整
  後補上）皆為常設回歸測試，任一階段的修改都不得破壞這四個極限。

### Edge / Abuse cases

| 情境 | 處理方式 |
|---|---|
| 新增的 G_i(X,F,Y) 在某個極限下退化成除以零（如 G5,X 在 X=0） | 模組 B 的鏈式法則自動化必須保留符號形式，不在推導階段代入具體數值，除零風險留到 H4 數值抽樣階段用有限差分/取極限方式處理，並在報告中明確標註哪些數值範圍不適用 |
| 兩條獨立推導路徑因為化簡順序不同，符號上不相等但數值抽樣相等（浮點誤差內） | 不得視為通過；H1/H2 要求符號等價性檢查優先於數值近似，數值抽樣只是輔助證據，不能取代符號比對（呼應你規格 H1 原文「不能只做內部自洽檢查」） |
| cadabra2 版本升級後，D 模組某條 workaround 不再必要，但另一條規則仍需要 | 逐條規則獨立管理狀態（見 `BugWorkaround` 實體），不得整批移除或整批保留 |
| 論文轉錄（`paper_equations/`）本身可能有誤植（如 `schutz_sorkin_vector.py` 曾發生過的式號誤植事件） | 轉錄依據必須是論文官方 LaTeX 原始碼或官方勘誤，禁止使用 OCR 或記憶轉錄；轉錄後應有獨立的人工核對步驟再進入驗收流程 |
| 協作 agent 在沒有 WSL 或沒有 cadabra2 的環境下嘗試執行 cadabra 相關階段 | Stage 計畫（見下）中每個涉及 cadabra 的階段，前置條件都要明確寫出「需要 WSL Ubuntu + cadabra2 2.5.14」，避免 agent 誤判環境已就緒 |

### 風險與未決事項

1. **【2026-09-18 部分解決，見 Stage 3】B2（δΓ 原語）在 L3（單層巢狀）上已確認可行。**
   `cadabra/variation_engine.py` 的 `delta_gamma_contracted_with()` 已在 cadabra2 2.5.14
   上機械化驗證：對 L3 = G3(X)∇_μA^μ 的 g^{μν} 變分，Palatini 展開、D1-safe 多項組合、
   δg_(下標)→δg^{上標} 轉換全部在 cadabra 裡正確跑通，最終閉式解
   `-1/2 G3,X A_μA_ν(∇_λA^λ)` 通過兩條獨立路徑驗證（G3=const 消失的 sanity check；
   以及對一個完全通用、非 FLRW 的度規做逐點有限微擾的獨立 SymPy 計算，符號差為 0——見
   `sympy_layer/delta_gamma_check.py`）。過程中發現兩個新的、與既有三個 bug 不同的 cadabra2
   摩擦點：(1) `canonicalise()`/`substitute()` 不會自動發現「交換一對啞指標」這種等價關係，
   必須手動驗證後當作已知恆等式使用；(2) `eliminate_metric()` 在純 `python3` 腳本（非
   TeXmacs/notebook 內核）下即使照搬官方 repo 自己的 worked example 也不會做出文件宣稱的化簡，
   這個坑目前繞過（改用明確的 `substitute()` 規則 + 手算/SymPy 交叉驗證完成最後的度量縮併與
   分部積分步驟），還沒有像既有三個 bug 一樣整理成 `bug_reports/` 最小重現案例。
   **L4（兩層巢狀）、L5（三層巢狀）是否還撐得住，仍然未知**——L3 只需要一個 δΓ，L4/L5
   需要巢狀 δΓ（δΓ 本身的變分又要再用一次 Palatini 展開），複雜度不是線性增加。建議下一次
   投入 L4/L5 之前，先評估這兩個新摩擦點（尤其 eliminate_metric 失效）在更深巢狀下會不會
   讓手動收尾的工作量爆炸；若真的撐不住，你規格 D 段提到的備案（改用 xAct/xTensor 產生獨立
   ground truth）需要認真考慮。
2. **【2026-09-18 更新，已在 Stage 2 解決】1703.09573 與 1603.05806 背景理論的差異範圍。**
   本規格原本假設兩篇論文的差異在於 G2 多了 F,Y 依賴、L5 多了 g5(X) 項，需要把這兩項加進
   背景推導。Stage 2 執行時直接下載了 arXiv:1703.09573v2 的官方 LaTeX 原始碼（不是 PDF
   OCR），發現論文本文明確寫著（astrophv2.tex 第 323-326 行）：「The functions g5, G6 and
   the additional dependence of F and Y in the function G2 ... do not contribute to the
   background equations of motion as expected.」也就是說：**g5、G6、以及 G2 對 F,Y 的額外
   依賴，對背景方程式（be1-be3，對應你原規格的 (2.11)-(2.13)）完全沒有貢獻**——不需要、
   也不應該被加進背景推導裡（F, Y 在這個純時間 ansatz 下本來就恆為零，這件事已經在
   `sympy_layer/proca_L6_background.py` 的檢查裡間接驗證過）。所以原本假設的「兩篇論文背景
   理論有差異」是不成立的：在背景層級，1603.05806（c2=0, d2 任意）與 1703.09573 的 L2-L5
   結構完全一致。真正需要修的不是「補項」，而是 `proca_minisuperspace.py` 本身三個獨立的
   既有 bug（EL_N 漏了 d/dt(∂L/∂Ṅ) 項、EL_a 漏了 Ostrogradsky 二階修正項、L5 的 G_munu
   nabla^mu A^nu 項用了一個實際上不成立的 Bianchi 恆等式化簡捷徑）——這三個 bug 都已在
   Stage 2 用 H1（逐字比對論文 be1-be3）與 H2（`proca_background_direct_route.py` 提供
   的獨立直接協變路徑，兩者完全符號相等）修復並雙重驗證，見
   `verification/registry.json` 的 `1703_09573_be1/be2/be3` 三筆記錄（皆為 `paper_matched`）。
   g5、G6、F、Y 的依賴仍然只在**微擾層**（純量/向量部門，Stage 4-5 範圍）才會真正出現，
   背景層不受影響。
3. **cadabra2 的三個 bug 未回報上游，長期維護風險。** 若未來 cadabra2 大改版（如 3.x），
   這三個 workaround 是否還適用是未知數。建議依本規格第 0 節的建議，把最小重現案例整理成
   issue 回報，即使不指望短期修復，也能在未來版本更新時有明確依據可查。
4. **階段五的「全新理論」示範案例尚未選定。** 建議延後到階段一至四完成後才選定，因為選錯
   案例（例如選了一個實際上需要模組 B 尚未支援的特殊變分技巧的理論）會讓階段五變成臨時抱佛腳。

---

## 非技術規格文件

> 這份文件是寫給不需要讀程式碼、但想知道「這套系統在做什麼、做到哪裡、能不能信任它算出來的
> 東西」的人看的——例如論文共同作者、審稿人溝通用的說明，或是你自己想要一份不用回頭翻程式碼
> 就能掌握全貌的摘要。

### 這個工具能做什麼

這是一套幫忙做「重力理論數學推導」的計算工具。物理學家研究修改重力理論時，常常需要從一個
很長的數學公式（作用量）出發，透過一步一步的數學運算（變分），推導出這個理論實際預測的
方程式。這個過程如果純手算，非常容易在某一步漏掉一項、或是搞錯一個正負號——而且錯誤往往
要等到很後面才會被發現，那時候前面所有基於錯誤結果的推導都要重做。

這套工具的目標，是讓這整個推導過程由程式自動完成，並且每一步都可以被檢查、被信任。你只需要
把想研究的理論的數學公式寫進去，工具就會算出這個理論實際的預測方程式，而不需要你自己一項
一項手推。

### 你會怎麼使用它

目前這是透過在電腦上執行一段一段程式來操作，不是點滑鼠操作的畫面工具。實際使用時的過程是：

1. 先把想研究的理論的數學式子，用工具規定的統一寫法寫下來（例如「這個理論多了一項跟磁場
   強度平方有關的耦合」）。
2. 執行推導。工具會把中間的每一步都印出來，讓你可以隨時停下來檢查某一步對不對。
3. 工具會用兩種不同的算法各推導一次同一個結果，並自動比對這兩次結果是否完全一致。如果不一致，
   工具不會自己猜哪個對，而是把兩邊的完整過程都攤開給你看，讓你自己判斷。
4. 如果你想拿這個理論已知的某一篇論文的公式來核對，工具可以把論文的公式輸入進去，跟自己算出
   來的結果做逐項比對，並且清楚告訴你「這是跟論文第幾條式子比對的」。

### 你會看到哪些主要畫面

因為這是命令列工具（沒有滑鼠可以點的視窗），你看到的是文字報告，格式大致長這樣：

```
[通過] 檢查 1：背景方程式的某一項應該完全化簡為零 ...... 通過
[通過] 檢查 2：兩條獨立算法的結果完全一致 ...... 通過
[尚未解決] 檢查 5：某個方程式理論上不該出現的項，目前還有殘留，原因待查 ...... 未通過
共 6/7 項檢查通過
```

重點設計是：**沒有解決的問題會被清楚標示出來、放在最顯眼的位置**，不會被其他已經通過的檢查
淹沒或掩蓋。這是為了避免「大部分檢查都通過了，就忽略那一小項還沒解決的問題」這種常見的疏忽。

### 操作流程

- 第一步：確認要研究的理論公式已經用統一的寫法登記好，避免不同人（或不同次使用）各自用不同
  的符號習慣，導致同一個物理量被定義成不一樣的東西（這是這套工具過去真的發生過的問題，曾經
  因為這樣導致一個正負號算錯，而且錯了一段時間才被抓出來）。
- 第二步：執行推導，中間結果會逐步印出。
- 第三步：如果要跟某篇論文的公式核對，先把論文的公式輸入進去（工具會要求你附上這是論文的
  第幾條式子，而且必須是照論文原文抄錄，不能是憑印象默寫或用不可靠的轉檔方式取得）。
- 第四步：工具比對推導結果跟論文公式是否完全相符。相符的話會被登記成「已驗證」；不相符的話
  不會被登記成已驗證，而是保留在「未解決」的清單裡等你進一步處理。

### 你會看到的提示語

| 情境 | 你會看到的文字 |
|---|---|
| 成功比對 | `[通過] 這條方程式已對照論文第 X.Y 條式子，並用兩條獨立算法互相驗證一致` |
| 兩條算法結果不一致 | `[不一致] 兩條算法得到不同結果，以下附上兩邊完整過程，需要人工判斷` |
| 卡住、暫時無解 | `[尚未解決] 卡在第 N 步，已嘗試過 A、B 兩種方法皆未成功，目前推測原因是 ...` |
| 抄自文獻、尚未獨立驗證 | `[未獨立驗證] 這個公式的形式來自文獻，尚未由本系統獨立推導確認` |

### 限制與注意事項

- 這套工具目前只能在特定電腦環境下執行（需要安裝特定的兩套數學計算軟體），不是隨開即用的
  網頁工具，也沒有計畫做成網頁工具。
- 目前這套工具還在建置階段：最核心的「自動幫任意理論做變分」的能力還沒有完全做出來，現在
  已經做出來的部分，大多是針對特定幾個理論、用比較直接（但一樣經過雙重驗證）的方式算出來的
  結果，還不是「你寫任何公式它都能自動處理」的完全體。
- 目前有一個已知但尚未解決的問題：某個方程式在加入某一項理論修正後，理論上應該完全化簡成
  一個簡單的代數式，但實際算出來還殘留了一個不該出現的項。這個問題已經被明確記錄下來，
  在還沒解決之前，不會被當作「已完成」對待。
- 有一部分物質（例如流體）的數學處理方式，目前用了一個計算上比較簡單、但嚴格來說不是原本
  要研究的那個物理模型的替代做法。這部分已知需要重做，才能真正符合最終要驗證的論文內容。

### 成功完成後會得到什麼

當這套工具完全建置完成後，你可以：直接寫下任何一個由重力場、額外物質場組成的理論公式，
系統就會自動算出這個理論實際預測的方程式，而且每一個結果都附有「怎麼算出來的」與
「已經用哪些方法驗證過」的完整紀錄，可以直接拿來核對或寫進論文，不需要你自己重新手算一遍
確認正確性。

### 常見問題與錯誤提示

- **「為什麼工具說兩條算法結果不一致，卻沒有直接告訴我哪個是對的？」** 因為工具設計上刻意
  不自己猜答案——如果兩條路徑算出不同結果，代表其中至少有一條有問題（甚至可能兩條都有問題），
  這需要人來判斷，工具只負責把兩邊過程完整攤開，不做自動選邊。
- **「為什麼有些結果明明數字對得上，卻還是被標成『尚未驗證』？」** 因為單純數字對得上（可能
  是巧合，或者只在某個特定情況下剛好對），跟「數學式子本質上完全等價」是兩件事。工具要求
  的是後者，這樣比較能避免「看起來對，其實只是碰巧」的誤判。

---

## Codex / Claude Code 分階段開發計畫

Stage 切分原則：以「可交付、可測試、可回滾」為單位，優先處理風險最高、最會影響後續全部工作的
項目（模組 A 地基、背景理論版本修正、B2 可行性驗證），G/F 模組的大重構排在確認 B2 可行之後，
避免在地基不穩時投入大量心力。以下每個 Stage 都同時提供 Codex 與 Claude Code 版本 instructions。

**【2026-09-18 執行進度】Stage 0-7 已依序自動執行完畢**（單一 Claude Code session，未中斷人工確認，
每個 Stage 完成後才進到下一個）。彙總結果：`py -3 verification/run_all_checks.py` 81/81 通過，
`py -3 -m pytest tests/` 91/91 通過，cadabra 層 4 支既有腳本（`eh_variation.py`/`fR_variation.py`/
`proca_L2.py`/`proca_L3.py`）與 3 支新增 bug 重現腳本在 WSL cadabra2 2.5.14 下全部正常執行。
過程中新增/修改的檔案清單、每個 Stage 的實際完成度與已知落差，見各 Stage 標題旁的日期註記，
以及本文件「驗收條件」章節的逐階段「目前狀態」欄。**這次執行過程中發現並修復了 4 個先前未被
發現的真實 bug**（詳見「風險與未決事項」與 Stage 2/3 段落）：`euler_lagrange_all()` 的 EL_N 漏了
d/dt(∂L/∂Ṅ) 項、EL_a 漏了 Ostrogradsky 二階修正項、L5 的 Bianchi 恆等式化簡捷徑本身不成立、
以及一個 φ=−A0 指標升降慣例混淆（這正是模組 A 存在的理由）——這些都是純內部自洽檢查
（原本 6/7 通過）從未能抓到、只有對照論文原文逐項比對（H1）才抓到的錯誤，直接驗證了你規格
「三、我對工作方式的要求」與模組 H1 的必要性。

### Stage 0：建立慣例地基與檢查骨架

- 目標：把模組 A（慣例單一來源）與模組 H 的一鍵檢查骨架建起來，讓後續所有 Stage 都有地基可站。
- 前置條件：無（可立刻開始）。

**Codex Instructions**
```text
[建議貼用方式] 直接貼給 Codex 執行；本專案沒有 AGENTS.md，若後續要新增，建議放在專案根目錄，
內容至少要包含「新增 cadabra 腳本前必讀 cadabra/cadabra_utils.py 的 D1-D3 規則」這條長期規則。

[任務範圍]
做：新增 conventions/gp_conventions.py（模組 A）；新增 verification/run_all_checks.py（彙總既有
12 支 sympy_layer 檔案與 4 支 cadabra 檔案的 run_checks()，此階段先支援 sympy_layer 部分，
cadabra 部分因需要 WSL 環境，先做成可選跳過）；新增 verification/registry.json 的初始 schema
（依 docs/SPEC.md 的「核心資料模型」章節的三個實體定義）。
不做：不修改任何既有理論推導邏輯；不改動 cadabra/ 既有四支檔案的推導內容（只允許之後 Stage
補齊 L3 的 g^{mu nu} 變分，這裡不動）。

[需修改/新增的檔案清單]
新增：conventions/__init__.py, conventions/gp_conventions.py, verification/__init__.py,
verification/run_all_checks.py, verification/registry_schema.py, tests/__init__.py

[具體步驟]
1. 讀 docs/SPEC.md 的「模組 A」「核心資料模型」「執行環境」三節。
2. 在 gp_conventions.py 定義：signature 慣例、t/k/Mpl 等共用符號、X/F/Y 的符號定義（不含具體
   ansatz 代入，只到符號層級）、L2..L6 的 G_i 函式符號登記。
3. 在 registry_schema.py 用 dataclass 定義 DerivedFormula / ConventionEntry / BugWorkaround
   三個實體（欄位依 docs/SPEC.md 核心資料模型章節）。
4. 在 run_all_checks.py 寫一個 runner：對 sympy_layer/ 下每支檔案的 run_checks() 呼叫並彙總
   結果，輸出文字報告（known_issue 排最前）與 JSON 摘要；cadabra/ 部分先跳過並在報告中標注
   "skipped: requires WSL"。
5. 執行 py -3 verification/run_all_checks.py，確認可以正確彙總既有 sympy_layer 全部檔案的結果
   （應該看到 proca_minisuperspace.py 的已知 open issue 被正確標記為 known_issue，不是被吞掉）。

[輸出格式要求] Python 檔案需含 module docstring 說明用途；registry.json 初始檔案為空陣列 []。

[測試要求] 用 py -3（Windows）執行 verification/run_all_checks.py，確認不崩潰且輸出包含至少
proca_minisuperspace.py 的已知 1 項 known_issue。

[驗收標準 DoD]
- gp_conventions.py 存在且可被 import。
- run_all_checks.py 能一次跑完全部既有 sympy_layer 檔案的 run_checks()，且已知的 open issue
  沒有被靜默吞掉。
- 沒有修改任何既有理論推導程式碼。
```

**Claude Code Instructions**
```text
[建議貼用方式] 直接貼給 Claude Code；建議把「新增 cadabra 腳本前必讀 D1-D3 規則」「sympy 用
py -3、cadabra 用 WSL」這兩條長期規則寫進 CLAUDE.md（若尚不存在則新建），供未來每個 session
自動載入，不必每次重新提醒。

[任務範圍] 同上 Codex 版本。
[需修改/新增的檔案清單] 同上 Codex 版本，另外新增/更新 CLAUDE.md（若已存在則追加，不覆蓋既有內容）。

[具體步驟]
1-5 同 Codex 版本。
6. 在 CLAUDE.md 追加一段「執行環境」說明：sympy 層用 `py -3`（Windows），不要用系統預設
   `python`/`python3`（沒有 sympy）；cadabra 層需要 WSL Ubuntu，指令前綴 `wsl.exe --`。
7. 在 CLAUDE.md 追加「新增 cadabra 腳本前，必讀 cadabra/cadabra_utils.py 開頭的 D1-D3 規則」。

[輸出格式要求] 同上。
[測試要求] 同上，另外確認 CLAUDE.md 的追加內容沒有破壞既有格式。
[驗收標準 DoD] 同上 Codex 版本，另外 CLAUDE.md 內含上述兩條規則。
```

- 風險與回滾方式：純新增檔案，無既有邏輯改動，回滾方式為直接刪除新增檔案。

---

### Stage 1：cadabra2 bug 上游回報 + D 模組擴大套用範圍

- 目標：把本規格第 0 節的調查結果落地成實際 issue 回報，並把 D 模組規則明確套用到後續所有新
  cadabra 腳本（用 lint 強制，不是靠自覺）。
- 前置條件：Stage 0 完成（有 run_all_checks.py 可以掛 lint 檢查）。

**Codex Instructions**
```text
[建議貼用方式] 直接貼給 Codex。

[任務範圍]
做：把 cadabra_utils.py 文件開頭描述的三個 bug，各整理成一個最小可重現的獨立 .cdb 或 .py 片段
（不依賴本專案其他程式碼），存到 cadabra/bug_reports/ 下三個檔案；在 verification/run_all_checks.py
中加入呼叫 cadabra/check_no_rename_dummies.py 的步驟，作為既有 lint 的一部分。
不做：不需要真的到 github.com/kpeeters/cadabra2 開 issue（那是人類決策，不是 agent 自主行為）；
只需要準備好可以直接貼上去的最小重現內容。

[需修改/新增的檔案清單]
新增：cadabra/bug_reports/accent_boundary_repro.py,
cadabra/bug_reports/rename_dummies_repro.py,
cadabra/bug_reports/substitute_dummy_collision_repro.py
修改：verification/run_all_checks.py（加入 check_no_rename_dummies 呼叫）

[具體步驟]
1. 讀 cadabra/cadabra_utils.py 全文 docstring，理解三個 bug 各自的觸發條件。
2. 為每個 bug 寫一個獨立、最小、不依賴 cadabra_utils.py 的重現腳本（因為要拿去給上游看，不能
   依賴本專案的 workaround 函式）。
3. 每個重現腳本要能在 WSL 的 cadabra2 2.5.14 下用 `wsl.exe -- python3 <script>` 執行並清楚印出
   「預期行為 vs 實際錯誤訊息」。
4. 更新 run_all_checks.py，在報告最前面加入 check_no_rename_dummies.py 的結果。

[輸出格式要求] 每個 bug_reports 檔案開頭要有：bug 一句話描述、預期行為、實際行為、cadabra2
版本。

[測試要求] 用 `wsl.exe -- python3 cadabra/bug_reports/<file>.py` 逐一執行，確認每個都能重現
docstring 描述的錯誤訊息。

[驗收標準 DoD]
- 三個最小重現腳本都能獨立執行且重現對應錯誤。
- run_all_checks.py 的報告最前面顯示 lint 結果。
- （人類步驟，不算入 agent DoD）你可自行決定是否要把這三個重現腳本貼到
  github.com/kpeeters/cadabra2/issues 開新 issue。
```

**Claude Code Instructions**
```text
[建議貼用方式] 直接貼給 Claude Code。
[任務範圍] 同上 Codex 版本。
[需修改/新增的檔案清單] 同上 Codex 版本。
[具體步驟] 同上 Codex 版本 1-4。
[輸出格式要求] 同上。
[測試要求] 同上，注意本機 bash 工具預設不在 WSL 內，需要明確用 wsl.exe -- 前綴執行。
[驗收標準 DoD] 同上 Codex 版本。
```

- 風險與回滾方式：純新增獨立重現腳本，無風險；回滾方式為刪除新增檔案。

---

### Stage 2：階段一核心——背景理論升級到 arXiv:1703.09573 真正版本

- 目標：把 `proca_minisuperspace.py` 從對照 1603.05806 升級成對照 1703.09573（G2(X,F,Y)、補
  g5(X) 項），完成你原規格的「階段一」驗收。
- 前置條件：Stage 0（有 registry 骨架可登記結果）；建議先完成「風險與未決事項」第 2 點（逐字
  比對兩篇論文第 II 節差異範圍）。

**Codex Instructions**
```text
[建議貼用方式] 直接貼給 Codex。

[任務範圍]
做：新建 verification/paper_equations/arxiv_1703_09573.py，把 (2.11)-(2.13) 逐字轉錄成 SymPy
表達式（依 docs/SPEC.md「H1：論文即標準答案」的要求，轉錄依據須為論文原文或官方勘誤，不可用
OCR 或記憶）；擴充 sympy_layer/proca_minisuperspace.py（或新建
sympy_layer/proca_minisuperspace_v2.py，若要保留舊版作對照），把 G2 從 G2(X) 擴充成 G2(X,F,Y)，
並在 L5 補上 -g5(X) F̃^{alpha mu} F̃^beta_mu nabla_alpha A_beta 項；跑 Euler-Lagrange，與
arxiv_1703_09573.py 的 ground truth 逐項比對。
不做：不需要處理張量/向量/純量微擾（階段二以後的範圍）；不需要重寫 L6（proca_L6_background.py
已經是對照 1703.09573 的正確版本，不需要改動，只需要在最終彙總時把它也納入同一套 registry）。

[需修改/新增的檔案清單]
新增：verification/paper_equations/arxiv_1703_09573.py,
sympy_layer/proca_minisuperspace_v2.py（若選擇不覆蓋舊檔）
可能修改：sympy_layer/proca_minisuperspace.py（若選擇直接升級原檔，需保留舊版本的 6/7 checks
邏輯作為 regression，不可刪除既有驗證過的部分）

[具體步驟]
1. 用 1703.09573v2.pdf 對照，確認 F, Y 的定義（F=-1/4 F_munu F^munu, Y=A^mu A^nu F_mu^al
   F_nu_al）與 g5(X) 項在 L5 中的確切位置（依你原規格提供的作用量定義）。
2. 把 (2.11)-(2.13) 轉錄進 arxiv_1703_09573.py，附上頁碼/式號、轉錄依據。
3. 擴充 build_lagrangian()：G2 改吃 X,F,Y 三個變數（F,Y 在純時間 A^mu=(A0(t),0,0,0) 背景下的
   具體值需要先推導——依 docs/SPEC.md 模組 C 的精神，這裡先用分量方式算出 F,Y 在此 ansatz 下
   的值，並標記這是過渡期做法，非長期目標的協變建構）；補上 g5(X) 項對背景 Euler-Lagrange
   方程的貢獻。
4. 重跑 run_checks()：既有 6/7 pass 的邏輯必須保留為 regression（不能因為新增項而破壞舊結果）；
   新增至少 3 項針對 F,Y,g5 項的獨立檢查。
5. 呼叫 H2 要求：用第二條獨立路徑（直接對協變場方程代入 ansatz，而非只走 mini-superspace
   Euler-Lagrange）交叉驗證至少一條新增的方程式，記錄進 registry.json。

[輸出格式要求] 每個新公式要標明 provenance；與論文比對的部分要標明對照到哪一條式號。

[測試要求] `py -3 sympy_layer/proca_minisuperspace_v2.py`（或升級後的原檔）執行 run_checks()，
既有 regression 全過，新增檢查至少覆蓋 G2 的 F,Y 依賴與 g5 項對背景方程的貢獻。

[驗收標準 DoD]
- verification/registry.json 中，對應 (2.11)-(2.13) 的三條 DerivedFormula 達到至少
  self_consistent 狀態（若能在此 Stage 內完成雙路徑驗證與論文比對，達到 paper_matched）。
- 既有 6/7 checks 的 regression 全部保留且通過。
- 不得刪除或靜默略過任何一項既有檢查。
```

**Claude Code Instructions**
```text
[建議貼用方式] 直接貼給 Claude Code；此 Stage 涉及較多物理判斷（F,Y 在特定 ansatz 下的值、
g5 項的正確位置），建議在動手前先用一句話跟你確認轉錄與擴充方向，而不是直接大改。

[任務範圍] 同上 Codex 版本。
[需修改/新增的檔案清單] 同上 Codex 版本。
[具體步驟] 同上 Codex 版本 1-5，額外要求：步驟 1 完成後，先停下來用簡短摘要跟你確認「F, Y,
g5 項的具體形式與位置」，取得你的確認後才繼續動手擴充程式碼（呼應 docs/SPEC.md「我對工作方式
的要求」第 2 點：卡住或不確定時先問，不要自行簡化物理設定）。
[輸出格式要求] 同上。
[測試要求] 同上。
[驗收標準 DoD] 同上 Codex 版本。
```

- 風險與回滾方式：若升級後既有 6/7 checks regression 失敗，優先回滾到 Stage 0 完成時的版本
  （git revert 到本 Stage 之前），不強行修補以求通過。

---

### Stage 3：B2 δΓ 原語可行性驗證（在 L3 的 g^{mu nu} 變分上）

- 目標：判斷協變引擎的核心（δΓ 原語）能不能撐住巢狀 ∇，用風險最低的案例（L3 只有一層巢狀）
  先驗證，避免直接挑戰 L4/L5 造成大量時間浪費。
- 前置條件：Stage 1（D 模組 lint 已可用）；需要 WSL 環境。

**Codex Instructions**
```text
[建議貼用方式] 直接貼給 Codex；此 Stage 需要 WSL 環境執行，Codex 若無法存取 WSL，本 Stage
應改由 Claude Code 版本執行，或由人類在 WSL 內手動跑最終驗證。

[任務範圍]
做：新增 cadabra/variation_engine.py，實作一個通用的 delta_gamma() 函式（依你原規格 B2 的定義：
delta(nabla_mu V_nu) = nabla_mu(delta V_nu) - delta_Gamma^rho_mu_nu V_rho，
delta_Gamma^rho_mu_nu = 1/2 g^rho_sigma (nabla_mu delta_g_sigma_nu + nabla_nu delta_g_sigma_mu
- nabla_sigma delta_g_mu_nu)）；套用在 proca_L3.py 目前明確留白的「L3 對 g^{mu nu} 的變分」上，
補完這個案例。
不做：不嘗試 L4/L5（留待下一個 Stage，取決於本 Stage 的可行性結果）；不修改 proca_L3.py 既有
已驗證的 A_mu 變分部分。

[需修改/新增的檔案清單]
新增：cadabra/variation_engine.py
修改：cadabra/proca_L3.py（在既有檔案基礎上，補上 g^{mu nu} 變分部分，並更新檔尾的「deferred」
說明為「已完成」或「仍未完成＋具體卡點」）

[具體步驟]
1. 讀 cadabra/proca_L3.py 檔尾的完整 NOTE（已經寫明 delta_Gamma 展開的具體需求）。
2. 在 variation_engine.py 實作 delta_gamma() 為可重複呼叫的通用函式，套用 D1-D3 規則（新
   dummy pair 一律用 fresh_indices()；含 \\delta{...} 的項一律逐項建構逐項替換）。
3. 用這個函式對 nabla_mu A_nu（L3 內出現的協變導數）套用一次 delta_gamma()，得到含
   delta_Gamma 項的完整展開。
4. 把結果代入 L3 = G3(X) nabla_mu A^mu 的 g^{mu nu} 變分，完成 proca_L3.py 原本留白的部分。
5. 記錄本次驗證花費的步驟數與是否有撞到任何超出既有 D 模組三個已知 bug 之外的新異常行為
   （若有新異常，比照 Stage 1 的模式整理成新的 bug_reports 條目）。

[輸出格式要求] variation_engine.py 的 delta_gamma() 函式必須有清楚的 docstring 標明對應你
原規格 B2 的數學定義。

[測試要求] `wsl.exe -- python3 cadabra/proca_L3.py` 執行，確認新增的 g^{mu nu} 變分部分能
完整跑完且輸出非平凡結果（不是報錯中止）。

[驗收標準 DoD]
- delta_gamma() 存在、可重複呼叫、有 docstring。
- proca_L3.py 的 g^{mu nu} 變分部分從「deferred」變成「完成」，或者明確記錄「嘗試過但卡在
  具體哪一步、原因是什麼」（不允許沉默失敗或跳過）。
- 若本 Stage 判斷 δΓ 原語在此案例上可行，在 Stage 4 前先向你回報「B2 可行性確認：建議下一步
  是否直接挑戰 L4」，讓你決定要不要提高本 Stage 之後的投入優先度。
```

**Claude Code Instructions**
```text
[建議貼用方式] 直接貼給 Claude Code；這是本規格判定風險最高的 Stage（B2 是否可行直接決定
後續策略），建議完成後主動總結「可行/不可行/部分可行＋具體卡點」回報給你，不要只留在程式碼
註解裡。

[任務範圍] 同上 Codex 版本。
[需修改/新增的檔案清單] 同上 Codex 版本。
[具體步驟] 同上 Codex 版本 1-5。
[輸出格式要求] 同上。
[測試要求] 同上，並在完成後用一段簡短文字（不超過 200 字）總結可行性判斷，供你決定後續投入
優先度（呼應 docs/SPEC.md「我對工作方式的要求」第 2 點）。
[驗收標準 DoD] 同上 Codex 版本。
```

- 風險與回滾方式：若 δΓ 在此案例上完全不可行，不強行湊出結果；回報具體卡點與已嘗試方法，
  依你原規格 D 段的要求，此時應該考慮 xAct/xTensor 備案，而非繼續在 cadabra2 裡硬湊——這個
  決策點必須交還給你，不是 agent 自主決定切換工具鏈。

---

### Stage 4：模組 G 重做——真正的 Schutz-Sorkin 純量流體

- 目標：用 `schutz_sorkin_vector.py` 的模式，替換掉 `scalar_plus_fluid_sector.py` 裡違反你
  規格 G 段要求的 k-essence 對偶做法。
- 前置條件：Stage 2 完成（背景理論已是正確版本，微擾層才有正確地基可站）。

**Codex Instructions**
```text
[建議貼用方式] 直接貼給 Codex。

[任務範圍]
做：新增 sympy_layer/schutz_sorkin_scalar.py，仿照 schutz_sorkin_vector.py 的模式，從
Schutz-Sorkin 作用量 S_M = -int d^4x [sqrt(-g) rho_M(n) + J^mu(d_mu ell + A_i d_mu B_i)] 出發，
針對純量型微擾（delta ell、delta J^0 等），推導出 delta rho_M 與速度勢 v 的定義與二階作用量
貢獻。
不做：不刪除 scalar_plus_fluid_sector.py（保留作為對照/歷史記錄，但在檔案開頭加註「此檔案的
k-essence 對偶做法不符合 docs/SPEC.md 模組 G 要求，已由 schutz_sorkin_scalar.py 取代，僅保留
供參考」）；不修改 schutz_sorkin_vector.py 本身（它已經是正確範本）。

[需修改/新增的檔案清單]
新增：sympy_layer/schutz_sorkin_scalar.py
修改：sympy_layer/scalar_plus_fluid_sector.py（僅加註淘汰說明，不刪除內容）

[具體步驟]
1. 讀 schutz_sorkin_vector.py 全文，理解它如何從 J^mu 出發推導 (3.20)-(3.23)。
2. 對純量部門，J^0 的擾動 delta J（對應純量密度擾動）與 ell 的擾動（速度勢）需要被納入——
   這是 schutz_sorkin_vector.py 背景設定裡明確設為零的部分（"delta J=0: no scalar density
   pert., v=0: no scalar velocity potential"），現在要反過來把它們打開。
3. 仿照 schutz_sorkin_vector.py 的 build_action_from_scratch() 模式：用顯式 eps 標記展開到
   二階，不要用 sp.series()（依開發重點第 8 條）。
4. 得到 delta rho_M 與 v 的定義後，與 scalar_plus_fluid_sector.py 現有的 (chiV, dphi_s) 二場
   系統做結構比對（不要求數值相同，因為底層物理模型不同，但要確認新的 (chiV, delta rho_M/v)
   系統同樣能給出非退化的 2x2 kinetic matrix）。
5. 更新 verification/registry.json，把純量物質部門相關的 DerivedFormula 從（若曾被誤登記為）
   paper_matched 降級為對應新檔案的正確狀態，並標明舊版 k-essence 結果的 provenance 為
   "superseded_non_conforming"。

[輸出格式要求] 檔案開頭 docstring 需明確標注「這是真正的 Schutz-Sorkin 推導，非 k-essence
對偶，取代 scalar_plus_fluid_sector.py 的物質部門定義」。

[測試要求] `py -3 sympy_layer/schutz_sorkin_scalar.py` 執行 run_checks()，至少驗證 delta rho_M、
v 的定義與 Schutz-Sorkin 作用量的直接展開一致（H2 的第一條路徑）。

[驗收標準 DoD]
- delta rho_M 與 v 的定義來自真正的 J^mu/ell/B_i 變分，不含 k-essence P(X)=X 對偶的痕跡。
- scalar_plus_fluid_sector.py 保留但已加註淘汰說明。
- registry.json 正確反映新舊兩者的 provenance 差異。
```

**Claude Code Instructions**
```text
[建議貼用方式] 直接貼給 Claude Code。
[任務範圍] 同上 Codex 版本。
[需修改/新增的檔案清單] 同上 Codex 版本。
[具體步驟] 同上 Codex 版本 1-5。
[輸出格式要求] 同上。
[測試要求] 同上。
[驗收標準 DoD] 同上 Codex 版本，另外要求：完成後用一段文字向你確認「新的 delta rho_M, v 定義
是否符合你預期的 Schutz-Sorkin 純量部門物理內容」，因為這是規格 G 段明確要求不能用替身的核心
項目，值得多一層人工確認（不要求你回覆才能繼續下一個 Stage，但要求明確留下待確認記錄）。
```

- 風險與回滾方式：若真正 Schutz-Sorkin 推導導出的 kinetic matrix 結構與現有 k-essence 版本
  差異過大導致無法在合理時間內完成，先保留兩版並存，在 registry.json 中都標記清楚 provenance，
  不強行讓其中一版消失；決定是否要投入更多時間完成，交還給你判斷。

---

### Stage 5：模組 F 收斂——統一的 SVT ansatz 建構器與通用 N×N kinetic matrix 萃取

- 目標：把張量/向量/純量三個微擾部門收斂到同一個 ansatz 建構器（F1/F2），並把 F6 從目前
  `scalar_plus_fluid_sector.py` 裡寫死的 2x2 特例，改成通用 N×N 萃取器。
- 前置條件：Stage 4 完成（有真正的 Schutz-Sorkin 純量部門可以餵給新的通用萃取器）。

**Codex Instructions**
```text
[建議貼用方式] 直接貼給 Codex。

[任務範圍]
做：新增 sympy_layer/svt_builder.py（F1：一次產生論文 (3.1)-(3.2) 全部微擾變數 alpha, chi,
V_i, h_ij, delta phi, chi_V, E_j，並提供 SVT 投影算子 F2）；新增 sympy_layer/kinetic_matrix.py
（F6：輸入任意 N 個場的二階作用量，輸出 N×N kinetic matrix 與 gradient matrix）；把
vector_sector.py、adm_scalar.py、scalar_sector.py、schutz_sorkin_scalar.py（Stage 4 產物）
改為呼叫 svt_builder.py 取得各自子集合的微擾變數，而非各自獨立手推。
不做：不重新推導既有已驗證的物理結果本身（只重構「怎麼取得微擾變數」與「怎麼萃取 kinetic
matrix」這兩個共用零件，既有的物理內容不變）。

[需修改/新增的檔案清單]
新增：sympy_layer/svt_builder.py, sympy_layer/kinetic_matrix.py
修改：sympy_layer/vector_sector.py, sympy_layer/adm_scalar.py, sympy_layer/scalar_sector.py,
sympy_layer/schutz_sorkin_scalar.py（改為 import 共用建構器，而非各自重複定義 abar/chibar 等
微擾變數符號）

[具體步驟]
1. 讀 vector_sector.py, adm_scalar.py, scalar_sector.py 三支檔案，列出它們各自定義的微擾變數
   符號（如 abar, chibar, dphi, chiV 等），確認命名不一致的地方（如 abar/chibar 對應
   alpha/chi）。
2. 在 svt_builder.py 用論文 (3.1)-(3.2) 的統一命名（alpha, chi, V_i, h_ij, delta_phi, chi_V,
   E_j）建構全部微擾變數，並提供依張量/向量/純量部門取子集合的介面。
3. 在 kinetic_matrix.py 實作通用 N 場版本：輸入一個二階 Lagrangian 與 N 個場的時間導數符號
   清單，回傳 N×N 的 K_ij = (1/2) d^2L/d(qdot_i)d(qdot_j)（既有 pert_engine.py 的
   kinetic_coefficient() 是 N=1 特例，此函式是其推廣，不是取代——pert_engine.py 保留）。
4. 逐一改寫既有三支檔案，把符號定義改成從 svt_builder.py import，並重跑各自的 run_checks()
   確認 regression 全過（symbol 改名不應改變任何數值/符號結果）。
5. 用 kinetic_matrix.py 重跑 Stage 4 產出的 (chiV, delta_rho_M 或 v) 二場系統，確認能得到與
   scalar_plus_fluid_sector.py 舊版結構類似（但物理內容正確）的 2x2 矩陣，並可平行擴充成未來
   N>2 的情況而不需要修改 kinetic_matrix.py 本身。

[輸出格式要求] svt_builder.py 的變數命名需與論文 (3.1)-(3.2) 完全一致，並在 docstring 附上
對照表（舊符號 -> 新統一符號）供之後查閱。

[測試要求] py -3 執行修改後的 vector_sector.py, adm_scalar.py, scalar_sector.py，確認既有
run_checks() 結果與修改前完全一致（純重構，不改變物理內容）；py -3 對 kinetic_matrix.py 新增
至少一個 N=3 的合成測試案例（非物理真實案例即可，只驗證函式本身的正確性）。

[驗收標準 DoD]
- 三個既有微擾部門檔案的 run_checks() 結果與重構前逐項相同（regression 零破壞）。
- kinetic_matrix.py 能處理 N>2 的合成案例，不是寫死的 2x2。
- svt_builder.py 是三個部門共用的唯一微擾變數來源。
```

**Claude Code Instructions**
```text
[建議貼用方式] 直接貼給 Claude Code；這是大範圍重構（修改三支既有已驗證檔案），建議先完成
svt_builder.py 與 kinetic_matrix.py 的獨立單元測試，確認新零件本身正確後，才動手改既有檔案，
降低一次性大改動出錯風險。

[任務範圍] 同上 Codex 版本。
[需修改/新增的檔案清單] 同上 Codex 版本。
[具體步驟] 同上 Codex 版本 1-5，額外要求：步驟 4 逐一改寫既有檔案時，每改完一支就立刻跑一次
該檔案的 run_checks() 確認 regression 沒破壞，再繼續改下一支，不要一次改完三支才測試（降低
出錯時難以定位是哪一支檔案造成的問題）。
[輸出格式要求] 同上。
[測試要求] 同上。
[驗收標準 DoD] 同上 Codex 版本。
```

- 風險與回滾方式：這是本計畫中對既有已驗證程式碼改動範圍最大的 Stage；若任何一支既有檔案的
  regression 在重構後失敗，優先個別回滾該檔案（保留其他已成功重構的檔案），不整批回滾，
  因為三支檔案的重構彼此獨立，沒有必要因為一支失敗而放棄全部。

---

### Stage 6（倒數第 2）：整合測試／回歸測試／邊界測試補齊 【2026-09-18 完成】

實際完成內容（88 項 pytest 測試全過，`py -3 -m pytest tests/ -v`）：
`tests/conftest.py`（sys.path 設定）、`tests/test_regression.py`（把
`verification/run_all_checks.py` 彙總的全部 sympy_layer `run_checks()`
結果，逐條轉成獨立 pytest 參數化案例，共 78 條；名稱含
"KNOWN ISSUE"/"KNOWN OPEN ISSUE" 的失敗會轉成 `xfail` 而不是直接
failure，若這類項目意外變成 xpass 或非預期失敗都會在 pytest 摘要裡清楚
顯示）、`tests/test_lint.py`（`test_no_rename_dummies_violation`）、
`verification/known_limits.py`＋`tests/test_known_limits.py`（H3 已涵蓋
的 7 個極限測試彙總索引，GR/GR+Maxwell/標準 Proca/f(R) 皆有，標準 Proca
背景極限與 Horndeski 對應極限標記為 `covered=False`——誠實記錄尚未做，
不是省略不提）、`tests/test_matter_sector_provenance.py`（
`test_schutz_sorkin_scalar_no_kessence_shortcut`：靜態掃描程式碼本體
（排除 docstring，避免把「解釋自己沒有用 k-essence」的說明文字誤判為
違規）確認沒有退回 k-essence 對偶；另外確認 `scalar_plus_fluid_sector.py`
仍保留 SUPERSEDED 標記）。cadabra 層本身的符號推導腳本（需要 WSL +
cadabra2，且是 print-based 展示腳本、不是 `run_checks()` 函式）未被
包進 pytest，維持原本 `wsl.exe -- python3 cadabra/<file>.py` 個別執行。



- 目標：把 `tests/`（目前為空）補成真正的回歸測試網，把「測試案例」章節列出的項目全部落地，
  並確保 Stage 0-5 累積的全部 regression 都被自動涵蓋，不再只靠人工個別執行 `run_checks()`。
- 前置條件：Stage 0-5 全部完成。

**Codex Instructions**
```text
[建議貼用方式] 直接貼給 Codex。

[任務範圍]
做：在 tests/ 下建立 pytest 測試，把 docs/SPEC.md「測試案例」章節列出的全部項目（
test_A0_algebraic, test_L6_background_zero, test_dual_route_consistency,
test_no_rename_dummies_violation, test_schutz_sorkin_scalar_no_kessence_shortcut,
test_known_limits）逐一實作為可獨立執行、可重複的 pytest 測試函式；補上 H3 的常設極限測試集
（verification/known_limits.py，收斂目前散落在各檔案裡的 GR/標準 Proca/f(R)/Maxwell 極限
測試）。
不做：不新增任何理論推導內容，本 Stage 純粹是把既有結果包成自動化回歸測試。

[需修改/新增的檔案清單]
新增：tests/test_background.py, tests/test_perturbation.py, tests/test_lint.py,
tests/test_known_limits.py, verification/known_limits.py

[具體步驟]
1. 逐一把「測試案例」章節的 6 個項目對應到既有程式碼的哪個函式/檔案，寫成 pytest wrapper。
2. known_limits.py 統一收斂目前散落在 proca_tensor.py（GR 極限）、adm_scalar.py（真空極限）、
   scalar_plus_fluid_sector.py（GR+fluid-only 極限）等檔案裡各自獨立的極限測試，改成單一
   可查詢的極限測試集合（不刪除原本檔案裡的極限測試，作為 regression 來源，known_limits.py
   是彙總入口）。
3. 用 `py -3 -m pytest tests/` 執行全部測試（cadabra 相關的測試如
   test_no_rename_dummies_violation 需要標記 `@pytest.mark.wsl`，並在 CI/本機說明如何用
   `wsl.exe -- python3 -m pytest tests/ -m wsl` 單獨執行 WSL 相關測試）。
4. 確認新的 pytest 測試網能夠正確地把 known_issue（如 A0 方程式殘留 addot）標記為
   `xfail`（預期失敗，而非直接刪除或忽略），避免這個已知問題在自動化測試裡被靜默通過或
   造成整個測試套件失敗兩種都不對的極端。

[輸出格式要求] 每個測試函式需要有 docstring 標明對應 docs/SPEC.md「測試案例」章節的哪一項。

[測試要求] `py -3 -m pytest tests/ -v` 全部通過或明確 xfail（已知問題），不允許有非預期的
失敗或錯誤。

[驗收標準 DoD]
- docs/SPEC.md「測試案例」章節列出的 6 項全部有對應的 pytest 測試。
- known_issue 用 xfail 明確標記，不是被跳過或刪除。
- pytest 測試網可以用單一指令執行（sympy 部分），WSL 相關部分有清楚標記的獨立執行方式。
```

**Claude Code Instructions**
```text
[建議貼用方式] 直接貼給 Claude Code；建議把「pytest 測試網的執行方式（sympy 用 py -3 -m
pytest，WSL 相關測試用 wsl.exe -- python3 -m pytest -m wsl）」補進 CLAUDE.md，讓未來每個
session 都知道怎麼跑測試，不需要重新摸索。

[任務範圍] 同上 Codex 版本。
[需修改/新增的檔案清單] 同上 Codex 版本，加上更新 CLAUDE.md。
[具體步驟] 同上 Codex 版本 1-4，額外在 CLAUDE.md 補上測試執行方式說明。
[輸出格式要求] 同上。
[測試要求] 同上。
[驗收標準 DoD] 同上 Codex 版本，另外 CLAUDE.md 已補上測試執行方式。
```

- 風險與回滾方式：純新增測試骨架，理論上不影響既有邏輯；若某個既有結果在包成 pytest 後才
  發現其實有問題（例如 regression 邏輯本身有誤），優先記錄為新的 known_issue，不強行讓測試
  通過而扭曲既有已驗證結果。

---

### Stage 7（最終）：可追溯性收尾與規格對齊文件化

- 目標：完成模組 I（LaTeX 輸出能力），並讓 `docs/SPEC.md` 的驗收條件章節與
  `verification/registry.json` 的實際狀態同步，交付一份「目前真實進度」的文件化成果。
- 前置條件：Stage 6 完成。

**Codex Instructions**
```text
[建議貼用方式] 直接貼給 Codex。

[任務範圍]
做：新增 traceability/latex_export.py，能把 registry.json 中 state >= self_consistent 的
DerivedFormula 輸出成附有來源標注的 LaTeX 片段；更新 docs/SPEC.md 的「驗收條件」章節表格，
把每個階段的實際完成狀態（依 registry.json 真實內容，不是憑印象）填入一欄「目前狀態」。
不做：不新增任何理論推導內容，本 Stage 是收尾與文件化。

[需修改/新增的檔案清單]
新增：traceability/__init__.py, traceability/latex_export.py
修改：docs/SPEC.md（驗收條件章節補上「目前狀態」欄）

[具體步驟]
1. 讀 verification/registry.json 目前內容，統計每個階段（對應你原規格的五個階段）各有幾條
   DerivedFormula 達到 paper_matched / regression_locked。
2. latex_export.py 對每條 DerivedFormula 輸出格式：LaTeX 公式本體 + 一行來源注解
   （% Source: arXiv:1703.09573 Eq (3.5), verified via <route_a> and <route_b>,
   last verified <date>）；provenance=literature_unverified 的項目輸出時要有明顯的
   \\textbf{[UNVERIFIED]} 標記，不能被排版掉。
3. 更新 docs/SPEC.md 驗收條件表格，如實填入「目前狀態」欄（可能是「部分完成」，不要為了
   好看而誇大）。

[輸出格式要求] latex_export.py 輸出的 .tex 片段要能直接貼進論文草稿而不需要手動調整格式。

[測試要求] 對 registry.json 裡至少一條已 paper_matched 的項目跑 latex_export.py，確認輸出
包含正確的來源標注格式。

[驗收標準 DoD]
- latex_export.py 能正確處理 provenance=literature_unverified 的視覺標記。
- docs/SPEC.md 驗收條件表格的「目前狀態」欄如實反映 registry.json 內容，不誇大也不隱藏
  known_issue。
```

**Claude Code Instructions**
```text
[建議貼用方式] 直接貼給 Claude Code。
[任務範圍] 同上 Codex 版本。
[需修改/新增的檔案清單] 同上 Codex 版本。
[具體步驟] 同上 Codex 版本 1-3。
[輸出格式要求] 同上。
[測試要求] 同上。
[驗收標準 DoD] 同上 Codex 版本，並在完成後用一段文字總結：目前五個驗收階段各自的真實完成度，
交給你做最終確認（不由 agent 自行宣稱「專案完成」）。
```

- 風險與回滾方式：純文件與輸出工具，無風險；若「目前狀態」欄位的統計與你自己認知的進度不符，
  以你的判斷為準，回滾方式為直接修正表格內容。

---

## 附錄：與你原規格逐條對照

本規格已涵蓋你原規格「二、系統必須具備的功能模組」A-J 全部十項、「三、我對工作方式的要求」
五點（分別落地為：H1/H2 的雙路徑與論文比對要求、Stage 說明中反覆出現的「卡住先問」指示、
provenance 標記機制、known_issue 狀態的誠實呈現、Stage 0 起就要求測試網先於大重構），以及
「四、驗收標準」五個階段（見「驗收條件」章節表格，逐一標注與你原描述的差異，主要差異是
新增了「必須是 1703.09573 真正版本而非 1603.05806」這個此前程式碼裡沒有意識到的落差）。
