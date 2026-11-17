# LevelSet 隐式曲面网格化：采样尺度、包围盒边界与交点精化对结果的影响分析

分析对象：`src/sdf.cpp`（Marching Tetrahedra / BCC 网格实现）、`test/sdf_test.cpp`、`test/context_test.cpp`、`src/edge_op.cpp`、`src/hashtable.h`、`src/impl.cpp`、`include/manifold/manifold.h`、`include/manifold/linalg.h`。下文所有结论均标注源码或测试证据；无法从代码推出保证的地方明确标为「无保证」。

---

## 1. 算法与数据流概述

`Manifold::LevelSet`（`src/sdf.cpp:473`）→ `Impl::CreateLevelSet`（`src/sdf.cpp:485`）：

1. **参数校验**（`sdf.cpp:497-501`）：包围盒非有限、min>max、`edgeLength` 非有限或 ≤0 → `MakeEmpty(Error::InvalidConstruction)`。`tolerance <= 0` 被替换为 `+inf`（`sdf.cpp:503-505`），即默认不做精化迭代。
2. **网格建立**（`sdf.cpp:509-525`）：`gridSize = trunc(dim/edgeLength + 1)`，`spacing = dim/(gridSize-1)`。采用体心立方（BCC）双网格：`w=1` 点位于整数格点，`w=0` 点偏移 `-0.5*spacing`（`Position()`，`sdf.cpp:130-132`）。
3. **体素采样**（`sdf.cpp:537-542`）：对整个网格（含外圈 1 层 padding）调用 `BoundedSDF` 一次性求值并存入 `voxels` 数组。
4. **近表面检测与顶点吸附**（`NearSurface`，`sdf.cpp:229-298`）：找出有符号变化的格点，靠近表面的格点被「吸附」到表面上（`FindSurface`），移动量每轴限制在 `kS=0.25` 倍 spacing 以内（`sdf.cpp:32,282`），以防四面体翻转。
5. **边交点计算**（`ComputeVerts`，`sdf.cpp:314-353`）：未被吸附的穿越边用 `FindSurface`（ITP 求根，`sdf.cpp:154-189`）产生顶点。
6. **四面体剖分建三角形**（`BuildTris`，`sdf.cpp:355-440`）：每个格点负责周围 6 个四面体，查表 `TetTri0/1` 生成三角形。
7. **后处理**（`sdf.cpp:646-656`）：`CreateHalfedges` → `CleanupTopology`（`SplitPinchedVerts` + `DedupeEdges`，`src/edge_op.cpp:47-56`）→ `RemoveUnreferencedVerts` → `CalculateBBox` → `SortGeometry` → 法线。顶点吸附可能产生「亲吻顶点」（even-manifold），注释明确说明：不做吸附「保证 2-manifold」，做了则必须靠 `CleanupTopology` 后处理恢复（`sdf.cpp:270-274`）。

**结果保证的基线**：输出恒为合法 2-manifold 封闭网格（`CleanupTopology` 内有 `DEBUG_ASSERT(IsManifold())`，`edge_op.cpp:49`），且多次运行结果逐位一致（`TEST(SDF, Determinism)`，`test/sdf_test.cpp:213-225`）。**但这是「采样符号场的忠实网格化」，不是「真实零等值面的忠实网格化」**——区别正是本报告的主题。

---

## 2. 采样尺度（edgeLength → spacing）的实际语义

### 2.1 spacing 只是「近似」≤ edgeLength

文档自述 edgeLength 是 "Approximate maximum edge length"（`sdf.cpp:459`）。由 `gridSize = trunc(dim/edgeLength + 1)`、`spacing = dim/(gridSize-1)`（`sdf.cpp:512-518`）可推出：当 `dim/edgeLength` 略小于整数+1 时 spacing 被放大，最坏情形（`dim/edgeLength → 2⁻`）`spacing → 2·edgeLength`。例如 `dim=2, edgeLength=1.9` → `gridSize=2` → `spacing=2.0 > edgeLength`。**结论：把 edgeLength 当作特征尺寸下限来用时必须留 2 倍余量**（见 §7 接受标准）。

### 2.2 可检测性的判据是「格点落入」，不是「特征与网格相交」

一个正值区域（实体）被看到的充要条件是：至少一个 BCC 格点的 `sdf - level > 0`；一个负值通道（空隙）被看到的充要条件是：至少一个格点落入其中。符号变化只沿网格边检测（`NearSurface` 的 7 条边，`sdf.cpp:244-268`）。因此：

- **薄层/小特征**（厚度 < 局部格点间距且无任何格点落入正值区）→ 完全不可见，静默消失。
- **窄通道/窄缝**（宽度 < 格点间距且无格点落入负值区）→ 通道被「焊死」，两侧实体连通——连通关系错误但无任何报错。

测试证据：

- `TEST(SDF, Resize)`（`sdf_test.cpp:180-195`）：`Layers` 在 z∈(−0.25,0.25) 有一个真实薄层，但网格（spacing=1）没有任何格点取到正值（z=0 在边界被钳为 ≤0，z=0.5 处函数值恰为 0，z=1 处为 −1），该层被整体丢弃——结果 `bounds.min.z ≈ 1.5` 而非 0，genus=−8（9 个分离组件）而非包含第 10 层。**这是「真实特征被静默遗漏」的直接测试证据。**
- `TEST(SDF, SphereShell)`（`sdf_test.cpp:68-79`）：壁厚 0.005 的球壳用 spacing≈0.01 采样，结果是 genus≈14235±1000 的「瑞士奶酪」——薄壁没有被整体丢弃，但被穿孔成上万个把手。**薄特征不一定消失，也可能以拓扑完全错误但逐点误差很小的形式出现**（每个顶点仍在真实表面附近，连通关系却全错）。

### 2.3 内存与内置守卫

- 体素数组：`maxIndex = 2^(1 + Σ gridPow_i)` 个 double，`gridPow_i = ceil(log2(gridSize_i + 3))`（`ComputeGridPow`，`sdf.cpp:108-116`；`sdf.cpp:525,534`）。每轴向上取整到 2 的幂，最坏放大约 2³=8 倍，即 **体素内存 ≈ 8·2^(1+Σceil(log2(n_i+3))) 字节**。
- 哈希表：`min(2·maxIndex, 10·√maxIndex)` 项（`sdf.cpp:545-556`），每项另有 `vertPos`/`vertKey` 各 7 倍预留（`sdf.cpp:557-558`）；表满时整体翻倍重跑（`sdf.cpp:567-585`）。
- 内置守卫只有两个：单轴 `gridSize ≤ 2^28`、索引总位数 ≤ 64（`sdf.cpp:513-524`），失败返回 `Error::ResultTooLarge`（测试 `TEST(SDF, InvalidGrid)`，`sdf_test.cpp:59-65`）。**这两个守卫防的是整数溢出，不是 OOM**：2^28/轴对应的体素内存远超任何实际内存。`Vec<double> voxels(maxIndex)` 分配失败只会以异常形式表现，**不会**转换为任何 `Error` 状态。因此有限内存预算必须由调用方在调用前自行核算（公式如上），不能依赖库报错。

---

## 3. 包围盒边界行为：钳制与封口

### 3.1 机制

`BoundedSDF`（`sdf.cpp:138-150`）对采样值做两层修改：

- 网格外（padding 层，`boundDist < 0`）：**不调用用户 sdf**，直接返回 0.0；
- 边界层（`boundDist == 0`）：调用 sdf 但返回 `min(d, 0)`，即**强制为非正（外部）**；
- 内部：原样返回 `sdf(p) - level`。

效果：边界层格点恒为「外部」，任何穿出包围盒的曲面都会在边界处产生符号变化，被 marching 剖分封上一个贴在包围盒面上的盖。`FindSurface` 在一端值恰为 0 时直接返回该端点（`sdf.cpp:156-159`），所以封口顶点精确落在包围盒面上；`Bound()`（`sdf.cpp:134-136`）再把所有顶点钳进 `[bounds.min, bounds.max]`。

### 3.2 后果

- **结果等价于「{x ∈ bounds : sdf(x) > level} 这个被盒子裁剪并封口的实体」**，而不是「真实 level set 与盒子的交集的精确表达」。封口面是沿网格的平面盖，顶点精确在盒面上。
- 测试证据：
  - `TEST(SDF, Bounds3)`（`sdf_test.cpp:136-153`）：半径 1.2 的球放进 ±1 的盒子，结果 genus=0、bounding box 精确等于 ±1——球被盒面平切封口，输出是合法实体但**不是球**。
  - `TEST(SDF, Bounds)`（`sdf_test.cpp:94-113`）：`CubeVoid`（盒外为正）在 ±2 盒中生成 genus=−1 的空腔实体，bbox 精确等于给定 bounds——说明「实体一直延伸到边界」时输出 bbox 恒等于输入 bounds，**bbox 触及 bounds 面即为发生过裁剪/封口的可观测信号**。
- **贴近包围盒的曲面被「封口」意味着**：库把盒子当作实体边界的一部分。若调用方本意是「盒子只是采样窗口、曲面在窗外还有延续」，则得到的封口实体在边界附近的几何与拓扑（连通性、genus）都是窗口 artifacts，不是真实曲面的性质。若曲面恰好完全在盒外（盒内无正值格点），结果为空——与「盒内本无实体」不可区分（见 §5）。
- 边界层仍会调用一次用户 sdf（`boundDist == 0` 分支），但盒外点永不调用（`boundDist < 0` 直接返回）——sdf 只需在闭盒内有定义。

---

## 4. 交点精化（FindSurface / tolerance / 顶点吸附）

### 4.1 默认（tolerance ≤ 0）：纯线性插值

`tolerance ≤ 0 → +inf`（`sdf.cpp:503-505`），ITP 循环条件 `frac > check` 立即为假（`sdf.cpp:168`），直接返回两端点的线性插值交点（`sdf.cpp:188`），**不做任何额外 sdf 求值**。此时顶点位置误差取决于 sdf 沿边的非线性程度，**无任何误差界**。

### 4.2 tolerance > 0：沿边的一维根精化

ITP 迭代（`sdf.cpp:154-189`）把符号变化区间收窄到 `2·tolerance/边长`，即**顶点沿该边方向距「采样函数的符号变化点」不超过 tolerance**。注意三层限定：

1. 这是**沿边的一维**界，不是到真实曲面的三维 Hausdorff 距离；
2. 仅当 sdf 是 1-Lipschitz 的真距离场时，它才近似等于「到真实曲面的几何距离 ≤ tolerance」；文档明确说 sdf 「不需要是真距离，甚至不需要连续」（`sdf.cpp:455-457`），对非距离场该界只关于「被采样函数沿边的根」，无几何意义；
3. 精化**只作用于已检测到的穿越边**——它不能找回任何未被采样的特征。tolerance 与 edgeLength 因此不可互相替代（§7）。

### 4.3 顶点吸附对闭合性与局部误差的影响

`NearSurface`（`sdf.cpp:244-296`）中，若某邻居采样值幅度 `> kD·|自身值|`（`kD = 1/kS − 1 = 3`，`sdf.cpp:32-34`），该格点被吸附到表面上，其所有穿越边共享一个顶点，避免短边和劣质三角形；但「对向顶点」超过 `kMaxOpposed=3` 时放弃吸附（`sdf.cpp:275`），退回普通边交点。吸附位移每轴 ≤ `0.25·spacing`（防四面体翻转，`sdf.cpp:282`），位置再被 `Bound()` 钳入盒内。

- 对**局部几何误差**：吸附把顶点真正放到表面上（配合 tolerance），比纯边插值更准；但 `|val|` 比值判据隐含假设 sdf 幅度近似距离——对只有符号的函数（§6），`|val|` 全为 1，吸附条件 `1 > 3·1` 恒假，自动退化为纯边交点，行为仍然安全。
- 对**闭合性/流形性**：吸附引入 even-manifold「亲吻顶点」，由 `CleanupTopology` 的 `SplitPinchedVerts` + `DedupeEdges` 修复（`edge_op.cpp:47-56`；注释 `sdf.cpp:270-274`）。**最终输出的 2-manifold 性是有保证的**（debug 构建有断言）。

---

## 5. 空结果、封口与采样遗漏分别意味着什么

| 现象 | 机制 | 含义 | 可观测性 |
|---|---|---|---|
| 空结果 + `Status()==NoError` | 盒内没有任何格边符号变化 | **两种可能不可区分**：(a) 合法空——盒内确实无 level set；(b) 采样遗漏——真实特征存在但无格点落入正值区（`Layers` 的 z=0 层即此情形，`sdf_test.cpp:180-195`） | 仅看输出无法区分；需要调用方几何先验（§7） |
| 空结果 + `Status()!=NoError` | 参数非法 / 超内置上限 / 被取消 | **拒绝结果**：`InvalidConstruction`（`sdf.cpp:497-501`）、`ResultTooLarge`（`sdf.cpp:513-524`）、`Cancelled`（`sdf.cpp:492-495, 570-573, 651-654`；测试 `context_test.cpp:1057-1093`） | 明确可见，必须检查 `Status()` |
| 曲面贴近/穿出包围盒被封口 | `BoundedSDF` 边界钳制（`sdf.cpp:145-149`） | 输出 = 实体 ∩ 盒子 + 盒面平盖；几何在盒面附近是窗口 artifact；`Bounds3` 证明 genus/bbox 都反映盒子而非真实曲面 | bbox 触及 bounds 面即为信号 |
| 窄通道被焊死 | 负值区无格点落入 | 连通关系错误（应分离的组件合并、孔洞消失），无报错 | 无直接信号；需先验最小通道宽 vs spacing |
| 薄层/薄壁 | 正值区无格点落入 → 整层消失；部分落入 → 穿孔 | `SphereShell`：genus≈14235（`sdf_test.cpp:76`）；`Layers`：整层丢失 | genus/体积与预期不符是唯一信号 |
| 部分组件被捕获、部分遗漏 | 各组件独立满足/不满足格点落入条件 | 非空但**不完整**；`Resize` 测试实际就是 10 层丢了 1 层 | 与「合法结果」不可区分，需先验组件数/体积 |

**关键结论：「局部误差小」与「拓扑正确」是独立的。** `SphereShell` 每个顶点都在真实表面 ±tolerance 附近（局部误差小），genus 却差了上万（拓扑全错）。逐点精度由 tolerance/吸附保证（在 §4 的限定下）；拓扑正确性**只**能由「spacing 小于最小特征厚度与最小通道宽度」这一采样条件保证，且该条件库无法自检——它不知道真实曲面在哪。

---

## 6. 只有符号值、不连续值、非有限值时的保证与失败可见性

### 6.1 只有符号（±1 式函数）

文档明确允许（`sdf.cpp:455-457`）。可提供的保证：

- 输出仍是合法、封闭、确定性的 2-manifold，且忠实于**采样符号场**（marching 剖分对符号模式是精确的）；
- 顶点吸附自动关闭（§4.3），退化为纯边交点，行为安全；
- **不能提供**：任何顶点位置精度界（tolerance 失去几何意义，`FindSurface` 找到的只是符号跳变点）；任何关于真实零等值面位置的推断（符号场不携带距离信息）。

测试证据：`TEST(SDF, SineSurface)`（`sdf_test.cpp:197-211`）用不连续的阶跃函数（±1）得到 genus=38、体积/面积与解析值吻合——符号场采样足够密时结果可以很好，但这是测试选择的采样密度带来的，不是库的保证。

### 6.2 不连续值

同样明确允许。`FindSurface` 把不连续点当作根（二分/ITP 收敛到符号跳变处）。测试证据：`Layers`（`sdf_test.cpp:31-36`）是逐点不连续函数；其测试还展示了一个细节：**采样值恰为 0 的格点被视为「在表面上」**——`FindSurface` 遇到端点值为 0 直接返回该点（`sdf.cpp:156-159`），`Resize` 测试中层表面精确落在半整数格点（z=1.5、18.5）即由此产生。

### 6.3 非有限值（NaN / ±Inf）——无保证，且失败基本不可见

这是保证链上最薄弱的一环，逐条给出源码证据：

- **符号判定**：`SameSide` 用 `>` 比较（`sdf.cpp:204-206`），`NaN > 0` 为假 → NaN 格点被当作「外部」。NaN 区域周围会产生**虚假的符号变化和虚假表面碎片**。
- **顶点位置**：`FindSurface` 对 NaN 输入，`d0/(d0-d1)` 为 NaN，默认路径直接返回 NaN 位置；`Bound()` 经 linalg 的 `min(a,b)=a<b?a:b`、`max(a,b)=a<b?b:a`（`include/manifold/linalg.h:437-450`）：`max(NaN, origin)=NaN`，`min(NaN, upper)=upper`——**NaN 坐标被静默钳成包围盒上角点**，输出全部有限。
- **状态可见性**：`CreateLevelSet` 全程**没有** `IsFinite` 检查（对比 `Warp` 路径在 `impl.cpp:540-541` 会置 `Error::NonFiniteVertex`）。结果 `Status()==NoError`，顶点却可能粘在盒角上。
- `ComputeVerts` 中的 `std::isfinite(neighbor.distance)` 回退（`sdf.cpp:333-336`）是处理哈希表未命中（默认 `GridVert.distance = NAN`，`sdf.cpp:198`），**不是**对用户 NaN 的防护。
- ±Inf：`+inf` 被视为内部；`inf/(inf-d)` 产生 NaN，退化为上一条。

**结论：sdf 返回非有限值时库不提供任何保证，且失败不会通过 `Status()` 可见。** 调用方必须自行消毒（sanitize）sdf 输出，或事后检查网格顶点有限性 / bbox 合理性。这属于「缺少保证的输出」：不是合法空，也不是拒绝结果，而是**状态为 NoError 的不可信结果**。

---

## 7. 有限内存预算下的结果接受标准（建议）

### 7.1 调用前（预算核算）

1. 由 `gridSize_i = trunc(dim_i/edgeLength + 1)` 计算体素内存上界 `8 · 2^(1 + Σ ceil(log2(gridSize_i+3)))` 字节，加上哈希表上界 `min(2·maxIndex, 10·√maxIndex)` 项 × 约 300 字节（`GridVert` + 7×(vec3+uint64) 预留，`sdf.cpp:545-558`），与预算比较。**库的 `ResultTooLarge` 不能替代这一步**（§2.3）。
2. 预算不足时的唯一正确调整是放大 `edgeLength` 或缩小 `bounds`（内存 ∝ 体积 / edgeLength³）；**不存在**「保持 edgeLength 而用别参数省内存」的开关。

### 7.2 调用后（接受/拒绝判定）

按顺序检查，任一失败即拒绝：

1. `Status() == NoError`——否则是**拒绝结果**（参数/超限/取消），不是几何答案。
2. 若结果为空（`IsEmpty()`）：只有在调用方先验允许「盒内可以没有实体」时才接受为**合法空**；否则按「采样遗漏」处理——加密网格重试或报告未知。**空 + NoError 本身不携带任何信息。**
3. `BoundingBox()` 与 `bounds` 比较：bbox 触及 bounds 面 ⇒ 发生了裁剪封口（§3）。仅当「裁剪+封口」正是预期语义时接受；否则扩大 bounds 重跑。
4. 拓扑/几何诊断量与先验比对：`Genus()`、`NumVert()`、`Volume()`、`SurfaceArea()`——现有测试全部用这些做验收（`sdf_test.cpp:76,105,144,169-171,186,206-208,253`）。`GetEpsilon()` 给出后续布尔运算的容差基准（`manifold.cpp:382` 注释）。
5. 顶点有限性抽查（sdf 可能非有限时，§6.3）——库不代查。

### 7.3 不能相互替代的取舍

| 参数 | 买不到的东西 | 原因 |
|---|---|---|
| `tolerance` ↓ | 任何拓扑正确性 | 只精化已检测穿越边上的顶点位置（§4.2）；未采样的特征对它不存在 |
| `edgeLength` ↓ | 盒外几何 | 采样域仍被 bounds 截断封口 |
| `bounds` 扩大 | 分辨率 | spacing ≈ edgeLength，与 dim 解耦（gridSize 同比放大，`sdf.cpp:512-518`）；只花内存，不买精度 |
| `bounds` 缩小 | 采样密度 | 同上；只改变裁剪窗口 |
| `level` 偏移 | 分辨率 | 只是换等值面（inset/outset，`sdf.cpp:462-463`；`Blobs` 用 level=0.5，`sdf_test.cpp:249`） |

**唯一**能买拓扑召回率的是更小的 spacing（内存、时间按立方增长）；**唯一**能便宜买局部精度的是 tolerance（每个顶点多几次 sdf 求值，不占网格内存）。两者不可互换。

### 7.4 调用方必须提供的几何先验 / 诊断信息

库自身完全不掌握以下信息，接受标准必须以调用方先验形式注入：

- **最小特征厚度 t_min 与最小通道宽度 w_min**：要求 `spacing ≤ min(t_min, w_min)`；鉴于 §2.1 的 2 倍放大，应设 `edgeLength ≤ min(t_min, w_min)/2`。这是拓扑正确的**唯一**充分条件，且库无法自检。
- **包围盒先验**：真实 level set 是否被 bounds 严格包含（建议各向留 ≥1 个 spacing 余量）；若否，必须明确接受「裁剪+封口」语义。
- **存在性先验**：盒内是否允许为空（用于解释空结果，§5）；可辅以「在已知内部点处 sdf 求值」这类廉价探针。
- **sdf 性质**：是否 1-Lipschitz 真距离场（决定 tolerance 的几何意义）；是否可能返回非有限值（决定是否必须消毒/事后检查）。
- **期望值域**：体积/面积/genus 的量级，用于 §7.2 第 4 步诊断。

---

## 8. 无法保证的部分汇总（明示）

1. **拓扑忠实性**：genus、组件数、连通关系与真实 level set 的一致性——无任何保证；仅当 spacing 小于所有特征尺度时由采样定理式的外部论证成立，库不验证。
2. **非空性的逆否命题**：真实 level set 与盒相交 ⇒ 结果非空——**不成立**（`Layers` 反例）。
3. **空结果的因果**：空 + NoError 无法区分「合法空」与「遗漏」。
4. **默认顶点精度**：tolerance ≤ 0 时顶点位置无任何误差界。
5. **tolerance 的三维几何意义**：仅对 1-Lipschitz 真距离场成立；对符号场/不连续场只是沿边根的位置界。
6. **非有限 sdf 值**：无保证、无状态报告，NaN 顶点被静默钳到盒角（§6.3）。
7. **内存安全**：除 2^28/轴与 64 位索引两个防溢出守卫外无任何内存保护；分配失败不走 `Error` 通道。
8. **封口面的几何质量**：盒面盖沿网格走样，不是真实曲面的延拓；其面积/法线只反映盒子。
9. **`CleanupTopology` 修复的语义**：亲吻顶点拆分保证 2-manifold，但在病态三角化下拆分方式可能改变局部连通关系（`edge_op.cpp:44-56` 注释）——流形性有保证，与真实曲面的对应关系没有。
