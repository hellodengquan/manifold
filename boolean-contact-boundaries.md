# Manifold 实体布尔中的接触关系与容差边界

> 分析对象：本仓库（`task/34`，commit `923fb9e`）的 Boolean3（3D）与 Boolean2（2D）流水线。
> 结论都标注了源码位置；凡是「实测」数据，均来自本机用本仓库源码构建的 `libmanifold`（Release，`MANIFOLD_PAR=OFF`，macOS arm64，默认 `Quality::GetRelativePrecision() == 1e-12`）得到的探针程序，探针代码不进入仓库。

---

## 0. 结论摘要

1. **接触关系的判定不在布尔内核里，而在「符号扰动 + 容差 tie」里**。Boolean3 用 `TieTol = min(epsP, epsQ)`（`src/boolean3.cpp:35-37`）定义「相等」，再用 `Moved`/`Shadows` 把每个 tie 按操作类型的固定方向打破（`src/boolean3.cpp:99-139,167-170`）。Boolean2 用 `EpsilonFromScale` 加顶点合并 / 关联预分裂（`src/boolean2_predicates.cpp:50-55`，`src/boolean2.cpp:300-467,760-831`）。两套 epsilon 都**由绝对坐标尺度推出**，不是由输入特征尺寸推出。
2. **三种操作对接触的解释不同，且这种不同是「设计上的符号扰动方向」，不是几何推导**。共面重合在 Add 下内部共享面消失、在 Subtract 下等于空、在 Intersect 下等于 P（完全重合时）；仅边接触/仅点接触在三种操作下都**不焊接**（Add 结果分成两个分量）。2D 完全同构：共享边在 Add 下并成一条 6 点轮廓，在角点接触下保持 2 条轮廓。
3. **拓扑流形保证不能推出几何接触关系不变**。`IsManifold()` 只检查半边配对一致（`src/properties.cpp:102-107`），`Is2Manifold()` 只额外排除重复有向边（`src/properties.cpp:113-128`），二者都不检查 ε-validity / 自交 / 接触是否被保留；自交检查是独立的、只查 `2*epsilon_` 带、且只在 debug 路径触发（`src/properties.cpp:138-140`，`src/csg_tree.cpp:196`）。`PairUp` 的注释自己就把两者分开了：「配对选择对流形性保证是任意的，但必须按几何有效排序，否则本算法退化为启发式」（`src/boolean_result.cpp:284-302`）。
4. **同一容差策略在缩放/平移下并不可交换**。`epsilon_ = relativePrecision * bBox.Scale()`，而 `Box::Scale()` 是**最大绝对坐标**（`src/properties.cpp:408`，`include/manifold/common.h:323-329`）。因此平移会放大 ε，非均匀缩放会改变「特征厚度 / ε」的比值，链式布尔会因 bbox 增长而逐步变粗；结果 bbox 又反过来驱动 `RemoveDegenerates` 消除**历史**特征（`src/boolean_result.cpp:933,955`）。实测：把一个 1e-9 宽的缺口和一个 1e6 远处的方块并起来，缺口会被**追溯性消解**。
5. **`IsEmpty()` 既不能证明「正确为空」，也不能证明「特征合理消解」**。两者都是 `Status() == NoError`（`src/impl.h:145`）；实测 `1e12` 平移下的 1 单位立方体自身布尔得到空且 `NoError`。唯一的显式错误只有 `ResultTooLarge` / `Cancelled` / 输入 `status_`（`src/boolean_result.cpp:738-747,781-785`）。
6. **仅凭面数或体积相等会漏掉**：接触被焊成 1 个分量 vs 保持 2 个分量（体积相同、面数不同）、共面接缝被合并 vs 因 `faceID` 而保留（体积相同、面数不同）、薄腔被消解 vs 从未存在（体积/面积/面数/分量数**全部相同**）、以及 ε 带两侧产生同一组数值。验收必须是**关系式 + 差分扰动 + 溯源审计**，不能是标量相等。

---

## 1. 判定点地图：接触与容差在流水线的哪里被决定

### 1.1 3D（Boolean3）

```
输入 Impl(P,Q)
  │  epsilon_ = relPrec * bBox.Scale()                  properties.cpp:388-409
  ▼
Boolean3::Boolean3
  │  tol = TieTol(P,Q) = min(epsP, epsQ)                boolean3.cpp:35-37,666
  │  bbox 早退：Box(P.min-2tol, P.max+2tol) ∩ Q.bBox    boolean3.cpp:667-674
  │  expandP_ = (op == Add)                             boolean3.cpp:653-661
  ▼
Intersect12 (×2, forward/reverse)                        boolean3.cpp:504-560
  │  Kernel02/11/12 用 Moved + Shadows 决定 tie          boolean3.cpp:167-296,298-447
  │  交叉点用 Intersect()/Interpolate()，tol 内吸附       shared.h:69-146
  ▼
Winding03 (×2)                                           boolean3.cpp:562-649
  │  连通分量 flood fill，再对分量代表点做 Kernel02      boolean3.cpp:576-637
  ▼
Boolean3::Result                                          boolean_result.cpp:696-983
  │  winding → inclusion（c1,c2,c3）                     boolean_result.cpp:708-710,789-802
  │  numVertR==0 → 返回空 Impl                          boolean_result.cpp:831-833
  │  SizeOutput：面保留当且仅当 Σ|inclusion| > 0          boolean_result.cpp:142-148,154-158
  │  AppendPartial/New/WholeEdges + PairUp                boolean_result.cpp:304-503,284-302
  │  Face2Tri → ReorderHalfedges                         boolean_result.cpp:936-943
  │  CalculateBBox（重新定 ε）                            boolean_result.cpp:933 → properties.cpp:388-409
  │  RemoveDegenerates(nPv+nQv) + RemoveUnreferencedVerts boolean_result.cpp:955-956
  ▼
输出 Impl（面片 triangulate、法向、triRef、ε）
```

关键点：**eps 在输入构造时确定一次，在结果 bbox 上重新确定一次**；`RemoveDegenerates` 用的是**结果**的 ε，因此会回头清理输入侧的特征。

### 1.2 2D（Boolean2）

```
Polygons a,b
  │  eps = InferEps(a,b) = EpsilonFromScale(bBox.Scale())   boolean2.cpp:277-283
  │  (k_budget=1000, kAlphaCoeff=12.37, kU=2^-53)            boolean2.h:33-34
  ▼
MergeVerts（并查集，距离 ≤ eps）                            boolean2.cpp:300-467
RemapAndCollapse（端点合并到同一点则丢边）                  boolean2.cpp:470-480
BuildIncidenceLists（顶点落在边内部 eps 内 → 预分裂）        boolean2.cpp:760-831
SweepWinding：
   排列扫描（发现交叉、7.6.2 块规则、竖向合并）              boolean2_sweep.cpp:528-582,160-191
   绕数扫描（IsInside 决定保留）                            boolean2_sweep.cpp:114-126,337-351
OutEdgesToPolygons（保留的有向边 → 正则化环）               boolean2.cpp:202-262
```

与 3D 的**结构性差异**：2D 输出是正则化的 `Polygons`，**无法表达孤立的 0 维或 1 维特征**（`docs/Boolean2.md:135-139`）。因此 2D 的边接触在输出中必然要么被并成环、要么消失，不存在「保留一条悬空边」的选项。3D 则用 `Decompose()` 得到两个分量来表达同一件事。

---

## 2. 三类接触关系的代码级判定

### 2.1 共面/面重合（coplanar, 2D 支撑集）

判定路径：两个面在同一平面上时，`Kernel02`/`Kernel11` 的点-面、边-面谓词全部落入 `|p−q| ≤ tol`，`Shadows` 走「按移动量决定」分支：

```cpp
// src/boolean3.cpp:167-170
inline bool Shadows(double p, double q, double movedP, double movedQ, double tol) {
  return std::abs(p - q) > tol ? p < q : movedP - p < movedQ - q;
}
```

移动方向由 `Moved` 固定：

```cpp
// src/boolean3.cpp:108-115
static Moved P(const Manifold::Impl& m, double tol, bool expandP) {
  const double step = tol / 4;
  return Moved(m, Displacements(m), (expandP ? 1.125 : -1.125) * step,
               step * vec3(1.0 / 61, 1.0 / 89, 1.0 / 113));  // 轴相关偏移，防止斜率为 1 时抵消
}
static Moved Q(const Manifold::Impl& m, double tol) {
  return Moved(m, Displacements(m), tol / 4, vec3(0.0));
}
```

含义：
- P 在 Add 时**向外**扩（`+1.125·tol/4`），在 Subtract/Intersect 时**向内**收（`−1.125·tol/4`）；Q 永远向外（`+tol/4`）。所以**共面时 P 的面在 Add 下是最外层**，在 Intersect/Subtract 下被 Q 压过。
- `Displacements` 优先让顶点沿「同时把相邻面各向外推等距」的方向移动（`src/boolean3.cpp:47-91`），保证面保持平行；这解释了为什么平台/台阶类输入（大量精确共面）§2.4 的 tie 会一致地倒向一侧。
- `Shadows` 的注释明确指出**不能**用 `movedP < movedQ` 直接比较（`src/boolean3.cpp:162-166`），否则 `Precision2` 与 `PointWithinToleranceOfCoincidentEdge` 会破。这就是「接触语义由两套值（原值与移动值）共同决定」的源码依据。

2D 侧对应的判定是 `PolySetAdd` 的「同向求和、反向抵消、零多重度湮灭」：

```cpp
// src/boolean2_sweep.cpp:143-153
void PolySetAdd(PolySet2& ps, vec2 a, vec2 b, int64_t m) {
  if (a == b || m == 0) return;          // 零长/零多重度丢弃（6.5）
  if (kLexLess(b, a)) { std::swap(a, b); m = -m; }
  auto [it, inserted] = ps.try_emplace({a, b}, m);
  if (!inserted && (it->second += m) == 0) ps.erase(it);
}
```

共享边（方向相反、大小相等）在 `PolySetAdd` 中精确抵消；这正是 2D 共面接触不会留下内部边的机制。

### 2.2 仅边接触

「仅边接触」在 3D 中的本质是：交集是 1 维集合（测度 0 的面积），`Intersect12` 可能记录到交叉点，但**没有任何面的 `inclusion` 能累积到非零**，于是：

- `SizeOutput` 的 `keepFace` 是 `x > 0 ? 1 : 0`（`src/boolean_result.cpp:142-148`），所有面 `sidesPerFacePQ` 为 0，结果 `numFaceR == 0`；
- 若体积/顶点也全为 0，`numVertR == 0` 直接返回空（`src/boolean_result.cpp:831-833`）。

在 Add 下，边接触的输入本身各面都保留，但**不允许共享一条边**——因为一条边被 >2 个面共享就不再是 2-流形。代码的做法是保留两份独立顶点（`DuplicateVerts` 按 `|inclusion|` 复制，`src/boolean_result.cpp:54-66,838-846`），于是输出是**两个各自闭合的分量**。

2D 侧：顶点合并只合并 `≤ eps` 的顶点（`src/boolean2.cpp:300-467`）、边分裂把落在边内部的顶点切成子边（`src/boolean2.cpp:760-831`）；共享边上的两个端点完全相同，于是环沿边界走通，成为单条 6 点轮廓。

### 2.3 仅点接触

与边接触同理，但更极端：唯一可能的共享集是单个顶点（3D）/单个点（2D）。在 **2-流形**约束下，一个顶点被两个实体共享会形成 pinch（`SplitPinchedVerts` 专门处理 pinch，`src/edge_op.cpp:1119-1222`），所以内核不会把点接触焊起来；输出保持两个分量。2D 中角点接触同理，输出 2 条轮廓。

### 2.4 薄壁接近容差尺度

三条独立的阈值线，都在 ε 量级，但**语义不同**：

| 机制 | 阈值 | 源码 | 后果 |
| --- | --- | --- | --- |
| 短边坍缩 | `|e|² > ε²` 才跳过 | `src/edge_op.cpp:88-96` | 边长 ≤ ε 的边被 `CollapseDegenerate` 消除 |
| 共线坍缩 | `Colinear(edge)` | `src/edge_op.cpp:133-176` | 与邻面法向偏差 ≤ ε 的边被折叠（减少共面细分） |
| 塌缩合并容差 | `tolerance = max(tolerance, ε)` | `src/edge_op.cpp:374-375` | `Simplify()` 不可能比默认 ε 更细 |
| 2D 顶点合并 | 距离 ≤ `eps` | `src/boolean2.cpp:300-467` | 薄边两端合并 → `RemapAndCollapse` 丢边 |
| 2D 关联预分裂 | 顶点到边线距离² ≤ `eps²·|e|²` | `src/boolean2.cpp:515-531` | 薄片被拉成 T 型并强制过点 |
| 三角化短边 | `|e|²·4 < ε²` | `src/polygon.cpp:355-357,644-645` | 耳切时提前裁掉「短到必然退化」的耳 |

源码中的关键「边界」写法是**平方比较**，因此判定阈值本身还有一个隐含的 `√` 敏感度：

```cpp
// src/edge_op.cpp:89-96
if (!halfedge_.Valid(edge) ||
    length2(vertPos_[halfedge_.End(edge)] -
            vertPos_[halfedge_.Start(edge)]) > epsilon_ * epsilon_)
  continue;
CollapseDegenerate(...);   // 边长 <= epsilon_ 的边
```

**实测（本机）**：边长/厚度与 ε 的比值就是分水岭，且分水岭非常陡（不是渐变）：

```
cube ε=1e-12（边长 2）
  内部薄腔 (1,1,t) 从 cube 减去:
    t = 1.0·ε → tri 12, area 24, vol 8, comps 1   ← 薄腔被完全消解
    t = 1.2·ε → tri 24, area 26, vol 8, comps 2   ← 薄腔保留
  薄壳 (cube 2 − 内立方体 2−2w):
    w = 1.0·ε → 空
    w = 1.2·ε → tri 24, area 48, vol 2.9e-11, comps 2
  两个 2×2×2 立方体面接触，gap 扫描（ε_c = 1e-12, TieTol=1e-12）:
    gap ≤ 0.999·ε → 合并（12 tri / 1 comp）
    gap ≥ 1.001·ε → 分离（24 tri / 2 comp）
```

分水岭在 `1×~2× TieTol` 之间，取决于配置：同尺度立方体面接触实测 `0.999ε`；混合尺度（`epsA=1e-6, epsB=5e-13`）实测 `≈2.0 × min(eps)`。源码中允许大于 `tol` 的机制是两处 `2*tol` 余量（`src/boolean3.cpp:521-530` 明确注释「Ties shift comparisons by up to tol... Twice tol strictly contains that」，同款见 `:39-43`、`:625,668`），它与 `Shadows` 的严格不等号叠加，使「能被发现为 tie 的」范围宽于 `tol` 本身。对球/圆柱等斜向面，实测分水岭更大（intersect 到 1.08 ε），这是固定轴投影插值（`GetAxisAlignedProjection`，`src/shared.h:161-180`）的离散误差，不是另一个阈值。

---

## 3. 三种布尔操作如何解释接触关系

### 3.1 winding → inclusion 的映射

```cpp
// src/boolean_result.cpp:708-710
const int c1 = op == OpType::Intersect ? 0 : 1;
const int c2 = op == OpType::Add ? 1 : 0;
const int c3 = op == OpType::Intersect ? 1 : -1;
// :789-802
i03 = c1 + c3 * w03;     // P 顶点按 Q 的绕数
i30 = c2 + c3 * w30;     // Q 顶点按 P 的绕数
```

| 操作 | P 顶点保留条件 | Q 顶点保留条件 | 法向 | `expandP` |
| --- | --- | --- | --- | --- |
| Add | `1 − w03 > 0`（在 Q 外） | `1 − w30 > 0`（在 P 外） | 原样 | true |
| Subtract | `1 − w03 > 0` | `−w30 < 0`（在 P 内 → 保留并翻面） | Q 翻转（`invertQ`，`boolean_result.cpp:787,175-182`） | false |
| Intersect | `w03 > 0`（在 Q 内） | `w30 > 0`（在 P 内） | 原样 | false |

`inclusion` 的绝对值还决定**顶点复制份数**（重合区域绕数 >1）：`DuplicateVerts` 复制 `|inclusion|` 份（`src/boolean_result.cpp:54-66,838-846`），`CountNewVerts` 同样（`:80-102`）。这就是为什么「加三次同一个立方体」仍然体积为 1 而不是 3。

### 3.2 逐操作、逐接触类的解释（3D，实测）

测试体：边长 2 的立方体（vol 8, area 24, ε≈1e-12）。`Q` 是 `P` 的平移副本。

| 接触关系 | `P + Q` | `P − Q` | `P ^ Q` |
| --- | --- | --- | --- |
| **完全重合**（Q=P） | 12 tri / 1 comp / vol 8 | **空** | 12 tri / vol 8（= P） |
| **面重合**（Q 沿 x 平移 2） | 12 tri（无溯源）/ 20 tri（`AsOriginal`）/ 1 comp / vol 16 / area 40 | 12 tri / vol 8（= P） | **空** |
| **部分共面**（Q 沿 x 平移 0.5） | 12 tri / vol 10 / area 28 | 12 tri / vol 2 / area 12 | 12 tri / vol 6 / area 20 |
| **仅边接触**（Q 平移 (2,2,0)） | 24 tri / **2 comp** / vol 16 / area 48 | 12 tri / vol 8 | **空** |
| **仅点接触**（Q 平移 (2,2,2)） | 24 tri / **2 comp** / vol 16 / area 48 | 12 tri / vol 8 | **空** |

解释与源码对应：

- **Add + 面重合 → 内部面消失，但接缝是否保留取决于来源标签**：添加时 P 向外、Q 向外（`Moved::P` 的 `+1.125·step`、`Moved::Q` 的 `+step`），位于内部的重合面 `inclusion` 为 0，被 `keepFace` 丢弃（`boolean_result.cpp:144`）。体积恒为 16、面积恒为 40（外部面积），**但三角形数在 12 与 20 之间变化**：两个 `Cube()`（无 original 关系）实测 12 tri；用 `Cube().AsOriginal()` 并把同一 original 的平移副本送出，实测 20 tri。`TEST(Boolean, FaceUnion)` 的期望正是 `{{12, 20}}`（`test/boolean_test.cpp:681-689`），因为它开头就调用了 `AsOriginal()`（`:682`）。原因：`Continuous()` 在 `meshID/faceID` 不同处禁止合并共线边（`edge_op.cpp:239-244`）。**所以「面接触 Add 后是 12 还是 20 tri」不是几何结论，而是溯源标签的结论。**
- **Add + 边/点接触 → 2 个分量**：测度 0 的交集无法在 2-流形里「焊」而不 pinch，因此保留两份顶点。`TEST(Boolean, EdgeUnion)`（`test/boolean_test.cpp:692-696`）与 `CornerUnion`（`:706-710`）期望 `{{8,12},{8,12}}`，即每分量 8 vert / 12 tri、共 2 个分量。这与「接触不改变各输入的闭合性」一致，也与「接触不产生共享拓扑」一致。
- **Subtract + 完全重合 → 空**：`P − P` 时 P 向内收、Q 向外扩，内收面全在 Q 内，`i03 = 1 − 1 = 0` 全丢；与 `TEST(Boolean, SelfSubtract)`（`test/boolean_test.cpp:474-480`）一致。
- **Intersect + 面/边/点接触 → 空**：`c1=0,c3=1`，需要 `w03>0`；但收缩/扩张后重合面被推开，`i03=0`，`numVertR==0` 早退（`boolean_result.cpp:831-833`）。实测完全重合的 `P^P` 反而返回 P（因为 Intersect 只收缩 P、扩张 Q，重合时仍判为重叠——源码自己也承认这一点：「Technically Intersection should contract inQ, but doing it this way makes Split faster」`src/boolean3.cpp:659-661`）。**这是三种操作在重合情形下唯一的不对称**。
- **PartialCoincidentFaceUnion**（`test/boolean_test.cpp:1389-1413`）专门覆盖「部分共面 + 面重合被删」：断言体积相加、面积减去两倍重叠面，即内部共享面必须消失。这从测试侧固定了「共面重合在 Add 下被消解」的语义。

### 3.3 2D 的解释（实测）

单位正方形 A、B：

| 接触 | `A + B` | `A − B` | `A ^ B` |
| --- | --- | --- | --- |
| 共享边 | 1 contour / 6 pts / area 2 | 1 contour / area 1 | **0 contour** |
| 角点接触 | **2 contours** / area 2 | 1 contour / area 1 | **0 contour** |
| 完全重合 | 1 contour / area 1 | **0 contour** | 1 contour / area 1 |
| 部分共面（错开 0.5） | 1 contour / 8 pts / area 1.5 | area 0.5 | area 0.5 |

与 3D 同构：**共享边被正则在环里合并（6 点），共享点不能合并（2 条轮廓）**。绕数规则见 `IsInside`：

```cpp
// src/boolean2_sweep.cpp:114-126
case WindRule::Add:       return w > 0;
case WindRule::Intersect: return w > 1;
```

Subtract 通过给 B 负多重度再用 Add（`src/boolean2.cpp:285-292`），所以它和 Add 共用同一条边界规则。

---

## 4. 「结果为空」与「薄特征被消解」的区别

### 4.1 二者的机制完全不同

**结果为空**：`NumTri() == 0`（`src/impl.h:145`）。列在下面的六条来源中，只有最后两条是**显式错误**；前四条都以 `NoError` 返回，其中第 4 条可能是静默丢体：

| 来源 | 代码位置 | `Status()` | 是否错误 |
| --- | --- | --- | --- |
| 输入为空 | `boolean_result.cpp:749-759` | 继承输入（可能是错误） | 取决于输入 |
| bbox 无重叠（含 2·tol 余量） | `boolean3.cpp:667-674` | `NoError` | 否 |
| 真实交集为测度 0（接触） | `boolean_result.cpp:831-833` | `NoError` | 否 |
| `MakeEmpty(NoError)`：bbox 非有限或被抽稀，注释写「Decimated out of existence」 | `properties.cpp:404-407` | `NoError` | **可能否，也可能是丢体** |
| `ResultTooLarge` | `boolean3.cpp:706-709` → `boolean_result.cpp:781-785` | `ResultTooLarge` | 是 |
| 取消 | `boolean_result.cpp:761-775` | `Cancelled` | 是（显式） |

注意 `MakeEmpty(Error::NoError)`——「被抽稀到不存在」被当成功处理（`src/properties.cpp:404-407`）。实测在 `T=1e12`、边长 1 的立方体上（ε=1），`a+b`、`a^b`、`a-b`、`a+a` **全部返回空且 `NoError`**。这就是「空」不能作为成功信号的最直接证据。

**薄特征被消解**：结果非空（或整体为空），但输入中某个特征在输出里没有任何表示。机制是第 2.4 节的短边/共线/合并阈值，以及结果 bbox 驱动的 `RemoveDegenerates`。

### 4.2 二者在输出上不可区分（关键）

实测：边长 2 的立方体减去内部薄腔 `(1,1,t)`。

```
t = 1.0·ε : tri=12  area=24  vol=8  comps=1   ← 与「原始立方体」逐字段相同
t = 1.2·ε : tri=24  area=26  vol=8  comps=2   ← 薄腔存在
```

也就是说：
- 若调用方输入的是一个**本该有腔**的体，输出是 `tri=12/area=24/vol=8/comps=1`；
- 若调用方输入的就是**实心立方体**，输出也是 `tri=12/area=24/vol=8/comps=1`；
- 若布尔**错误地**丢掉了腔，输出还是同一组数。

**输出中没有任何残留痕迹**：没有 ε 之外的顶点、没有零面积面、`NumDegenerateTris()==0`、`Status()==NoError`。所以「合理消解」与「错误丢失」在单次调用的输出上**信息论上不可区分**，只能靠输入侧的先验（特征尺度 vs ε）或差分实验来判定。

### 4.3 更严重的一种：整体被消解

薄壳 `w = 1.0·ε` 的结果是**空**（vol 0, area 0），而 `w = 1.2·ε` 是 `comps=2`。一个原本体积为正的实体整体消失，`Status()` 仍是 `NoError`。同样地，独立存在的薄板（100×100×1e-11，自身 ε=5e-11）：

```
standalone:            tri=12, vol=1e-7, deg=0, 存在
plate + plate(副本):    空
plate.RemoveDegenerates(): 空
plate + 远处方块(无交):   comps=2, 板仍在（因为走 Compose 而非布尔）
plate - 角部探针:        空
```

**同一个网格「是否存在」取决于是否被一个布尔操作碰过**——这不是几何事实，而是流程事实。

---

## 5. 有效闭合体保留了哪些输入关系

### 5.1 保留的（源码可证）

1. **闭合性与定向**：所有有向边成对；`KeepFace` 只保留 `Σ|inclusion| > 0` 的面（`boolean_result.cpp:142-148`），`PairUp` 把起点/终点配对成半边（`boolean_result.cpp:284-302`），`Face2Tri` 三角化并 `ReorderHalfedges`（`:936-943`）。Subtract 时 Q 的面法向与 `backSide` 一起翻转（`:175-182,536-539`）。
2. **逐面来源关系 `triRef`**（meshID / faceID / triID）：`UpdateReference` 把输出的 `triRef` 映射回 P 或 Q 的原始 `triRef`，并给 Q 的 meshID 加偏移（`boolean_result.cpp:505-541`）。这是调用方唯一能用来做「哪些输入面活着」审计的东西。
3. **顶点属性**：`CreateProperties` 用重心坐标插值，并沿边界/顶点做 key 合并（`boolean_result.cpp:574-690`）；`GetBarycentric` 在 `ε` 内把点吸附到顶点/边（`src/shared.h:182-229`），所以属性在 ε 带内有量化。
4. **分量的拓扑量**：`Decompose()` 按有向边连通性重新分组（`src/constructors.cpp:433-492`），每分量重算 bbox 与 ε；`Genus()` 由 `1 − χ/2` 得到（`src/manifold.cpp:423-426`）。

### 5.2 **不**保留的（容易误以为保留）

1. **`OriginalID`**：布尔结果走完整路径时 `originalID` 是 `-1`（`UpdateReference` 只搬 `meshIDtransform`，不设 `originalID`；`InitializeOriginal` 只在构造/`AsOriginal`/hull/sdf 时调用，见 `src/impl.cpp:188`、`src/manifold.cpp:469`）。实测：`AsOriginal(11) + AsOriginal(22)` → `OriginalID() == -1`。但**空操作数短路径**会直接返回另一边（`boolean_result.cpp:749-759`），且结果**取决于空体的来源与是否已被求值**：

   | 空体状态 | `a + empty` | `empty + a` |
   | --- | --- | --- |
   | 默认构造、未求值 | 11 | **−1** |
   | 默认构造、已 `IsEmpty()` | 11 | **−1** |
   | 由布尔产生（`cube ^ 远处 cube`）、未求值 | 11 | 11 |
   | 同上、已 `IsEmpty()` | 11 | **−1** |

   几何（`vol/tri`）在四种情形下完全相同，**只有元数据不同**。根因是默认 `Box()` 的 `min=+inf, max=−inf`（`include/manifold/common.h:292-307`）让 `BatchUnion` 的「互不相交就 Compose」分组（`src/csg_tree.cpp:508-521`）随求值时机改变。**结论：`OriginalID` 不能作为溯源的可靠依据，必须用 `faceID`/`triRef`。**
2. **接触关系本身**：不存在「这两个分量在某个面上接触」的编码。`Decompose()` 之后你只看到两个独立分量；`MinGap()` 对接触与 1e-13 重叠都返回 0（`properties.cpp:477-495`，实测 `contact=0`、`overlap 1e-13=0`、`gap 1e-13≈1e-13`），因此 `MinGap` **能证明分离，不能证明接触**。
3. **共面接缝是否合并由 `faceID` 决定，而不是由几何决定**：`TriID2FaceID` 文档说「by making unique faceIDs, it prevents collapsing coplanar edges」（`src/manifold.cpp:514-528`）。`Continuous()` 用 `meshID`/`faceID`/属性索引判定两侧是否连续（`src/edge_op.cpp:239-244`），`Collinear()` 只在 `Continuous` 为真时才允许折叠共线边（`edge_op.cpp:141-146`）。实测：两个面接触立方体，默认并集 12 tri；给两侧不同 `faceID` 后并集 20 tri（接缝保留）。所以**同一条几何接缝，是否出现在输出里，取决于输入的来源标签**。
4. **`RemoveDegenerates(int firstNewVert)` 的参数**：文档（`src/edge_op.cpp:72-74`）声称「index < firstNewVert 的顶点跳过昂贵检查」，但函数体（`:79-131`）**完全没有使用该参数**，`grep firstNewVert src/` 只命中声明与注释。`boolean_result.cpp:955` 传入的 `nPv+nQv` 因此是**没有效果的**。这意味着重复布尔**不会**跳过旧顶点，历史特征会被反复重新评估并侵蚀——这正是 §7.2-L4 实测到的追溯性消解。

---

## 6. 拓扑流形保证是否足以推出几何接触关系不变？

**不足够。** 五条独立理由，均有源码依据：

1. **`IsManifold` 只查组合一致性**：半边数 %3、配对有效、`startVert` 合法（`src/properties.cpp:102-107`，`CheckHalfedges`）。`Is2Manifold` 只额外要求「不存在重复的有向边」（`:113-128`）。两者都**不看坐标**，因此对接触关系完全无话可说。
2. **自交检查是独立的、且不在生产路径**：`IsSelfIntersecting` 明确注释「not checking for epsilon-validity」，且只在 `selfIntersectionChecks` 开启时调用（`src/properties.cpp:134-140`，`src/csg_tree.cpp:196`），默认关闭。
3. **`PairUp` 自己把两者分开**：注释写「配对选择对 manifoldness 保证是任意的，但必须按几何有效排序；否则本算法退化为启发式」（`src/boolean_result.cpp:284-302`）。这等于承认：**流形性是组合不变量，几何有效性是另一个独立要求**。

4. **`RemoveDegenerates` 可以改变接触关系而保持流形性**：文档明说边坍缩可能「删除 handle 或分离网格，从而降低 Genus」，只有「塌缩成一对三角形后被整体移除」时 Genus 才增加（`src/edge_op.cpp:58-78, 1119-1222`）。也就是说 **genus / 分量数可以在保持流形的前提下变化**。
5. **符号扰动按操作类型分叉**：`expandP_ = (op == Add)`（`src/boolean3.cpp:655`），`Moved::P` 的符号随 `expandP` 变（`:108-112`）。同一个几何接触在 Add 与 Intersect 下被决定成相反结果，而两者输出都通过同样的 `IsManifold` 检查。

**实测反例（最短路径）**：

```
T=1e9, ε=1e-3, 边长 1 立方体:
   a + b(平移 1)      → tri=12, vol≈2, comps=1   （面接触被焊）
   a ^ b(平移 1)      → 空                         （接触判为无交）
   a + a              → tri=12, vol≈1              （自并被去重）
T=1e12, ε=1:
   a + b(平移 1)      → 空, Status=NoError
   a + a              → 空, Status=NoError
   1e12 处的 a 本身   → vol=0.99994（顶点已被量化）
```

同一个几何（边长 1、间距 1、面接触），在 T=0 时是「并集体积 2 的连通体」，在 T=1e12 时是「空」。两者都满足流形性（空网格 `IsManifold()==true`，`properties.cpp:103`；`IsManifold`/`Is2Manifold` 只在 `src/impl.h:161-162` 的内部 `Impl` 上，不是公开 API）。

更直接的对照（实测，间隙 `1e-4` 的面接触并集）：

```
T=0,    ε=5e-13, gap=1e-4  → comps=2, tri=24   （分离：gap > ε）
T=1e9,  ε=1e-3,  gap=1e-4  → comps=1, tri=12   （焊接：gap < ε）
```

同一对输入、同一个间隙，仅因整体平移就被判定成相反的接触关系。

**保证边界的准确表述**：

> Manifold 保证：输出是组合意义下 2-流形、定向一致的网格；每个输出三角形可追溯到一个输入 `triRef`；体积/面积按 Kahan 求和（`properties.cpp:300-330`）。
>
> Manifold **不**保证：输出与输入在 ε 带内的几何关系（接触/重叠/间隙）在任何意义下保持不变；也不保证同一操作在不同平移、缩放、链式顺序下给出相同的接触分类。

---

## 7. 输入缩放后采用同一容差策略的限制

### 7.1 ε 的定义是「绝对坐标尺度」，不是「特征尺度」

```cpp
// src/properties.cpp:408
epsilon_ = Quality::GetRelativePrecision() * bBox_.Scale();

// include/manifold/common.h:323-329
constexpr double Scale() const {
  vec3 absMax = la::max(la::abs(min), la::max(max));
  return la::max(absMax.x, la::max(absMax.y, absMax.z));
}
```

2D 同理，`InferEps` 用 `bBox.Scale()`，并**明确声明不是平移不变的**：

```cpp
// src/boolean2.cpp:272-276
// ... Position-inclusive (bBox.Scale() = max abs coordinate),
// ... so the arrangement stays resolvable far from origin
// ... without a recenter. Not translation invariant by design.
```

### 7.2 由此产生的五条限制

**(L1) 平移放大 ε，且是线性的。** 实测（边长 1 的立方体）：

| 平移 T | ε | ε/尺寸 |
| --- | --- | --- |
| 0 | 5e-13 | 5e-13 |
| 1 | 1.5e-12 | 1.5e-12 |
| 1e3 | 1.0e-9 | 1.0e-9 |
| 1e6 | 1.0e-6 | 1.0e-6 |
| 1e9 | 1.0e-3 | 1.0e-3 |
| 1e12 | 1.0 | 1.0 |

因此「同一几何」在原点附近和在 1e9 处会走不同的分支。实测（§6 的对照）：间隙 `1e-4` 的面接触并集，`T=0`（ε=5e-13）判为分离（2 comps），`T=1e9`（ε=1e-3）判为焊接（1 comp）。**同一组输入位移，结果随平移翻转。**

更强的一条：**同一世界坐标几何，因为"先切缺口"还是"先平移"的顺序不同，结果不同**（实测，1e-9 宽缺口）：

```
在原点切缺口，再平移到 1e6    → ε=1e-6, tri=28, area=6.08   （缺口存活：Translate 是惰性的，只重算 ε，不调用 RemoveDegenerates）
直接在 1e6 处切缺口           → ε=1e-6, tri=12, area=6.00   （缺口在建布尔时已按 ε=1e-6 消解）
直接切完再平移回原点           → ε=5e-13, tri=12, area=6.00  （已丢失，平移回来无法恢复）
```

这解释了为什么 **事后归一化无法挽回**——损失发生在构造期的布尔里。

**(L2) 均匀缩放在原点附近近似协变，非均匀缩放不协变。** 因为 `Scale()` 取三轴最大绝对值，压缩某一轴不降低 ε。实测：

```
100×100×1e-9 薄板:        ε=5e-11,  t/ε = 20
Scale(1, 1, 0.01) 之后:   ε=5e-11,  t/ε = 0.2   ← 特征变成亚 ε
Scale(0.01,0.01,1) 之后:  ε=5e-13,  t/ε = 2000  ← 特征变得很安全
```

同一块板，只因为先做了不同的缩放，从「安全」变成「必被消解」。

**(L3) `TieTol` 取两者较小者，与「各自的表示精度」不匹配。**

```cpp
// src/boolean3.cpp:35-37
inline double TieTol(const Manifold::Impl& a, const Manifold::Impl& b) {
  return std::min(a.epsilon_, b.epsilon_);
}
```

`a.epsilon_` 是 a **自身几何可信到**的量级。取 min 意味着：把一个在 1e6 处建模（ε=1e-6）的粗网格和一个在原点建模（ε=5e-13）的细网格做布尔时，tie 带是 5e-13，**远小于粗网格顶点的表示误差**。此时粗网格的「面」在 5e-13 尺度上并不可信，判定结果由舍入噪声主导。这不是「更精确」，而是**精度错配**。

实测（大立方体边长 1e6 于 x∈[−1e6,0] 与单位立方体相邻）：

```
eps(big)=1e-6, eps(small)=5e-13, TieTol=min=5e-13
gap 的合并/分离翻转点 ≈ 1.0e-12 = 2.0 × TieTol = 1e-6 × eps(big)
输出 ε（结果 bbox 决定）= 1e-6
```

也就是说：**判定用的容差是细网格的 5e-13，而结果的 `GetEpsilon()` 是粗网格的 1e-6**。调用方如果按「结果 ε 的 1e-6」判断「这两块应该被合并」，就会对 `1e-12` 与 `1e-6` 之间的一切情形做出错误预测。这就是精度错配的可观测后果：**tie 的尺度和结果的尺度不是同一个 ε**。

**(L4) 结果 bbox 重新定 ε，并驱动 `RemoveDegenerates` 追溯清理历史特征。** 这是最危险的一条：

```cpp
// src/boolean_result.cpp:933   （结果 bbox → 新 ε）
outR.CalculateBBox();
// :955-956
outR.RemoveDegenerates(nPv + nQv);
outR.RemoveUnreferencedVerts();
```

实测：一个带 1e-9 宽缺口、ε=5e-13 的立方体（缺口在自身 ε 下是真空的，`notchW=0`）。下表用 `NumTri` 与 `SurfaceArea` 作为判据——**不要**用 `WindingNumber` 判断合并后的网格，因为合并后 bbox 的 ε 也会污染查询（见 §R5；实测合并网格在同一点给出 `notchW=2`，这是查询侧 ε 所致，不是几何变化）：

```
fine                                 ε=5e-13  tri=28  area=6.0800   （缺口存在）
fine + 1e6 远处方块（无交，走 Compose） ε=1e-6   tri=40  area=12.0800  （几何不变）
  Decompose 后的 near 分量            ε=5e-13  tri=28  area=6.0800   （缺口仍在）
merged.RemoveDegenerates()            ε=1e-6   tri=24  area=12.0000  ← near 分量塌成 12 tri / 6.000
merged.Simplify()                     ε=1e-6   tri=24  area=12.0000  ← 同上
merged + 远处 tip（一次真正的布尔）     ε=1e-6   tri=40  area=12.0200  ← near 分量塌成 12 tri / 6.000
```

（`area` 是整网格所有分量之和；near 分量自身从 6.080 掉到 6.000，正是 1e-9 缺口消失造成的 0.08 表面积损失。）

即：**一个完全无关、相距 1e6 的布尔，会因为结果 bbox 横跨 1e6 而把 ε 抬高 6 个数量级，然后回头把近处的历史特征吃掉**。注意触发条件是**结果** bbox 覆盖近处——若布尔的结果只落在远处（或只落在近处），ε 不会被抬高，特征就安全（实测见 §9 C4 的三行对照）。`RemoveDegenerates(firstNewVert)` 本应避免这一点，但该参数未实现（§5.2 第 4 条）。

**(L5) 相对精度下限是 1e-12，无法收紧。**

```cpp
// src/manifold.cpp:57,69-71
static std::atomic<double> relativePrecision_ = kPrecision;   // kPrecision = 1e-12
void Quality::SetRelativePrecision(double p) {
  relativePrecision_ = std::max(p, kPrecision);               // 只能放松，不能收紧
}
```

`utils.h:40` 定义 `kPrecision = 1e-12`。所以调用方**无法**通过 `SetRelativePrecision` 表示「我的坐标是双精度精确的，请用 1e-15」。唯一办法是**把坐标归一化**（平移+缩放到单位尺度），做完布尔再变回来；但变回来时的取整又会引入新的、可能大于原 ε 的位移。

### 7.3 文档与实现不一致（应视为契约风险）

`docs/Boolean2.md:149-152` 写：「Inputs are translated into a local frame before the arrangement is built, then translated back on output.」但实现明确没有这样做：

```cpp
// src/boolean2.cpp:174-177
// No local-origin recenter: the inferred eps is position-inclusive
// (bBox.Scale()), matching Boolean3's MaxEpsilon and the polygon
// triangulator, so the arrangement stays resolvable far from origin.
```

实测 2D：共享边并集在 T=0 / 1e6 时面积=2（正确），在 **T=1e9 时面积=2.5**（错误：多出一个角点，6 点轮廓变成 5 点/6 点的不一致环）。这说明文档承诺的「局部坐标系」并不存在，且大坐标下的误差已经超出「舍入」而进入了「错误拓扑」。

---

## 8. 调用方如何区分「合理的特征消解」与「错误地丢失实体」

因为没有 API 直接回答这个问题，需要一套**输入侧先验 + 输出侧审计 + 差分实验**的组合。按可靠性从高到低：

### R1. 先归一化，再做布尔（最强的预防手段）

```
c = (bbox.min + bbox.max)/2
s = max(bbox.Size())                     // 或 max|coord|
P' = (P - c)/s ;  Q' = (Q - c)/s
R' = op(P', Q')
R  = s*R' + c
```
归一化后 `ε' = relPrec · O(1) ≈ 1e-12`。于是「特征相对厚度 ≥ ~1e-11」一定安全。

**必须在每个布尔之前归一化，而不是之后**。实测：

```
(world 坐标 1e6 处，1e-9 缺口)
  在原点建 + 切缺口 → Translate(1e6)         tri=28, area=6.08   （活）
  在 1e6 处直接建 + 切缺口                    tri=12, area=6.00   （已死）
  在 1e6 处建 + 切缺口，再 Translate(-1e6)     tri=12, area=6.00   （不可逆）
```

`Translate` 是惰性的（`src/manifold.cpp:537-539`），只在求值时经 `Impl::Transform` 重算 bbox/ε（`src/impl.cpp:549-617`），**不调用 `RemoveDegenerates`**——所以「先切后移」保住缺口，而「先移后切」在切的时候就按大 ε 消解了。

**代价**：`s*R' + c` 的往返变换有取整，可能把顶点移动 > 原 ε；只应在 `|c|/s` 远大于 1 且特征厚度也远大于 `relPrec*|c|` 时使用。若特征本身已经相对 `|c|` 小于 `relPrec`，归一化也救不回来——此时必须先降 `|c|`（改写模型坐标系）。

### R2. 用输入自身尺度做前置判断

对每个操作数算
```
feat = 最小正特征尺度  （min 非零边长、最小面间距、最小非共面顶点距离）
eps  = GetEpsilon()
```
- `feat > 4·eps` → 该特征**应该**存活；若输出中缺失，判为**错误丢失**。
- `feat < eps/4` → 消解在契约内，判为**合理**。
- `eps/4 ≤ feat ≤ 4·eps` → **歧义区**，任何结果都不可断言；调用方应拒绝或标注。

这个 `4×` 系数来自 §2.4 实测分水岭（≈1.0×）加安全边际；源码的比较都是 `> ε²` 或 `> tol` 的严格不等号，所以临界值本身不稳定。

### R3. 溯源审计（唯一的「谁活着」证据）

```cpp
MeshGL64 pm = P.GetMeshGL64(); pm.faceID.assign(pm.NumTri(), 100);
MeshGL64 qm = Q.GetMeshGL64(); qm.faceID.assign(qm.NumTri(), 200);
auto out = (Manifold(pm) - Manifold(qm)).GetMeshGL64();
// 期望：Q 的 200 只应出现在「被切开的新内表面」上；
//       若 Q 体量 >> ε³ 而 200 完全缺席 → Q 被整体丢失
```
实测（薄腔，`faceID=99` 是切削体，`7` 是被切的本体）：
```
t = 4.0·ε → 输出 faceIDs = {7, 99}   （切削面保留）
t = 1.2·ε → 输出 faceIDs = {7, 99}   （仍在）
t = 0.9·ε → 输出 faceIDs = {7}       （切削体整体消失，无错误）
```
这个阶梯是**可观测的**，可以在调用方建立「最小可表征厚度」的自检。注意 `OriginalID` 不可用于此（§5.2 第 1 条），只能用 `faceID`/`triRef`。

**实现提示**：`faceID` 是输入 `MeshGL` 的字段，`UpdateReference` 会把它搬到输出（`boolean_result.cpp:505-529`）。但 `faceID` 同时会影响共面边合并（`Continuous()`，`edge_op.cpp:239-244`），所以**给每个面一个唯一 `faceID` 会改变几何输出**（实测面接触并集从 12 tri 变 20 tri）。做溯源审计时应接受这一点，或在审计副本上与生产副本分开运行。

### R4. 差分扰动（检测「你是否站在 ε 边缘」）

对可疑特征，把某个参数（间隙、厚度、平移）在其标称值上做 ±k·ε 扫描（k=1,2,4,8）：

- 结果在该带内**保持同一组合类型** → 稳定。
- 在带内**翻转一次**且翻转点 ≈ `1~2·ε` → 你正好在 ε 边缘，应视为不可判定。
- 在带内**翻转多次**或出现非单调 → 实现缺陷或接触歧义。

实测的面接触并集就在 `0.999ε → 1.001ε` 之间翻转一次（同尺度）；混合尺度下翻转点 ≈ `2·min(eps)`（§7.2-L3）。两次都是**单次、可预期**的翻转。

### R5. 点成员判定（验证几何而非拓扑）

用 `WindingNumber()` 在**距任何面 ≥ 4·ε'** 的采样点上验证（`ε'` 是**该结果网格**的 `GetEpsilon()`），并与解析真值比较。注意：

```cpp
// src/boolean3.cpp:747,757,777
const Box tolBox(bBox_.min - epsilon_, bBox_.max + epsilon_);
pointImpl.epsilon_ = epsilon_;             // 查询方与网格同 ε
auto f = [&](int i) { return VerticalQueryBox(..., epsilon_); };
```

查询侧也用**该网格的 ε**。所以对大 bbox 网格做细特征查询，结果本身会被 ε 抹平。实测：带缺口的网格在自身 ε=5e-13 下查询得 `notchW=0`（缺口存在）；把它和一个 1e6 远处的方块组合后（合成网格 ε=1e-6），同一坐标查询得 `notchW=2`（缺口两侧壁被当作重合，各贡献一次穿越）。**因此 R5 必须在「目标特征所属分量的 ε」下做（先 `Decompose()` 再查），而不是在整结果 bbox 的 ε 下做。**

### R6. 明确区分 `Status()` 的语义

| 判据 | 能说明什么 | 不能说明什么 |
| --- | --- | --- |
| `Status() == NoError` | 没有数值/规模/取消错误 | 不代表任何几何正确性 |
| `IsEmpty()` | 结果为 0 三角形 | 不能区分「正确空」与「丢体」 |
| `NumDegenerateTris() == 0` | 无共线三角形 | 不代表无亚 ε 特征 |
| `MatchesTriNormals()` | 面法向与顶点顺序一致（容差 `1000·ε`，`properties.cpp:198`） | 不检查接触 |
| `Genus()` | `1 − χ/2`，仅单分量有意义（`manifold.cpp:423-426`） | 分量数变化不报错、且不是公开的流形性检查 |
| `MinGap() > 0` | **确实分离** | `== 0` 不能区分接触与亚 ε 重叠 |

---

## 9. 相容的验收判据

设 `ε = TieTol(P,Q) = min(epsP, epsQ)`，`d` 为两个实体的特征间距（重叠为负）。以下判据互相兼容（可同时通过），并且**任何一条都不能被单纯的标量相等替代**。

### A 组：核心不变量（必须始终成立）

- **A1 状态完整**：任何几何异常**不得**以 `NoError` + 空结果的形式静默出现。可判定的必要检查：
  - `op == Add`、两个操作数都非空 → 结果**不可能**为空。返回空即错误（A1-a）。
  - `op == Subtract`、两个操作数都非空 → 结果**不可能**比 P 大，但仍然可能因 Q 完全包含 P 而为空；此时须用 `WindingNumber` 抽样确认 P 是否真被完整覆盖，不能仅凭 bbox。
  - `op == Intersect` → bbox 重叠**不能**推出交集非空（非凸体可交错），因此不能仅凭 bbox 报错；应改用一组 `WindingNumber` 采样点（同时在 P、Q 内且距面 >4ε）作为「必须非空」的证据。
  - 当前实现在 `T=1e12` 下 `a+a`、`a+b` 皆空且 `NoError`（§4.1）；这是 A1-a 的直接违反。
- **A2 闭合性**：每个非空分量满足 2-流形 + 定向一致；`NumDegenerateTris()==0`；`Volume()>0`、`SurfaceArea()>0`；`Genus()` 为整数。2D 对应「每个环自闭合、边平衡（出度+入度）为 0」。**注意**：`Is2Manifold()`/`IsManifold()` 只在 `Impl` 上，**不是公开 API**（`src/properties.cpp:102-128`，`include/manifold/manifold.h` 仅暴露 `MatchesTriNormals/NumDegenerateTris/Genus/GetEpsilon`）；生产构建里 `boolean_result.cpp:945-965` 的 `DEBUG_ASSERT` 也被编译掉（`optional_assert.h:67-68`）。因此调用方若要 A2，必须自己从 `GetMeshGL()` 重建半边并检查，或自行链接内部头。
- **A3 能量守恒（弱）**：`|op 后的体积 − 解析体积| ≤ C·ε·Area`；2D 用面积。**只作必要非充分条件**。

### B 组：ε 带外稳定性（真正的接触判据）

- **B1 组合稳定性**：对 `|d| ≥ 4·ε` 的所有采样，输出的**分量数、每分量 genus、每分量体积、以及存活 `faceID` 集合**必须恒定。
- **B2 单调过渡**：在 `|d| ≤ 4·ε` 内允许结果变化，但**只允许一次**翻转（例如 2 comp → 1 comp），且不得出现分量的「出现—消失—再出现」。
- **B3 操作语义一致性**（对接触集 `d = 0`，区分「有面积重叠」「零测度接触」「完全重合」）：
  - **有面积共面重叠**（面部分重合，Q 相对 P 平移 `0 < t < 尺寸`）：Add → `comps == 1`、`Vol = Va + Vb − Vov`（内部共享面消失，`Area = Aa + Ab − 2·Aov`）；Subtract → `comps == 1`、`Vol = Va − Vov`；Intersect → `Vol = Vov`（**不是** `min(Va,Vb)`）。
    实测：边长 2 立方体沿 x 平移 0.5 → Add `vol 10 / area 28`、Subtract `vol 2 / area 12`、Intersect `vol 6 / area 20`（`Vov = 1.5·2·2 = 6`）。
  - **零测度面接触**（两面仅贴合）：Add → `comps == 1`、`Vol = Va + Vb`、`Area = Aa + Ab − 2·A_shared`；Subtract → `Vol = Va`；Intersect → `IsEmpty()`。
  - **零测度边/点接触**：Add → `comps == 2`、`Vol = Va + Vb`；Subtract → `Vol = Va`；Intersect → `IsEmpty()`。
  - **完全包含**（Q ⊂ P，无面接触）：Intersect → `Vol = Vb = min(Va,Vb)`；Subtract → `Vol = Va − Vb`。
  - **完全重合**（Q ≡ P）：Add → `Vol = Va`；Subtract → `IsEmpty()`；Intersect → `Vol = Va`。
    注意这条与「零测度面接触的 Intersect → 空」**不同**，是三种操作里唯一一处不对称，源码自述为 Split 的性能取舍（`src/boolean3.cpp:659-661`）。
  - **这条判据对应源码 §3.1-3.2 的 `c1/c2/c3` 与 `expandP`，是唯一能验证「三种操作对接触的解释是否正确」的判据。**

### C 组：缩放/平移协变性（针对 §7）

- **C1 均匀缩放协变**：对 `s > 1`，`op(sP, sQ) ≈ s·op(P,Q)`，容差 `4·ε·s`。
- **C2 平移协变（受限）**：仅当平移后 ε 不增大超过 4 倍时成立：
  ```
  要求  Scale(T(x)) ≤ 4·Scale(x)  即  |T| ≤ 3·|x|
  ```
  超出该范围，**必须**先归一化（R1），否则不要求协变，而应要求**显式报警**。
- **C3 非均匀缩放**：必须按「缩放后的 `feat/ε` 比值」重新评估，不得沿用缩放前的结论（实测 `Scale(1,1,0.01)` 把 `t/ε` 从 20 降到 0.2）。
- **C4 链式无侵蚀**：在结果中加入一个与近处分量距离 `> 4·ε` 的远处实体，**不得**改变近处分量的 `NumTri/Area/faceID 集合`。触发条件是「该布尔的结果 bbox 同时覆盖两者」，此时结果 ε 由大 bbox 决定，`RemoveDegenerates` 会回头清理近处（§7.2-L4）。实测对照，说明判据应以**结果 bbox** 为准：
  ```
  merged(=near + far@1e6) + tip@1e6（与 far 相交，结果跨 1e6）  → near 28 tri/6.08 → 12 tri/6.00
  merged + 新方块@2e6（与两者都不交，走 Compose）              → near 保持 28 tri/6.08
  merged - far（结果只含 near，bbox 回到近处）                  → near 保持 28 tri/6.08
  ```
  当前实现违反第一条，因此这条应作为**回归目标**，而不是现状断言。

### D 组：溯源完整性

- **D1 覆盖**：输出的 `faceID` 集合 ⊆ 输入的 `faceID` 集合（允许 meshID 偏移，见 `UpdateReference`）。
- **D2 必要面保留**：任一输入面，若其面积 `> 16·ε²` 且位于最终边界上，其 `faceID` **必须**出现在输出中。
- **D3 关系保留**：`OriginalID` 在走完整布尔路径时恒为 -1；仅在空操作数早退（`boolean_result.cpp:749-759`）或未求值 CSG 节点直接合并时可能保住输入值，且实测随空体的来源与求值时机变化（§5.2 第 1 条）。验收应**只断言「不依赖 `OriginalID`」**，把所有溯源检查放在 `faceID`/`triRef` 上。

### E 组：消解标注

- **E1**：若某输入特征 `feat < ε/4` 且输出中缺失，写入 `dissolved` 标注（当前 API 没有，需要调用方自行在 `faceID` 层面记录）。
- **E2**：若 `ε/4 ≤ feat ≤ 4·ε`，调用方**不得**对结果做任何几何断言；应视为「容差歧义」，类似于数值计算中的 `NaN`。

### 9.1 仅凭面数或体积相等遗漏了什么

| 被遗漏的维度 | 反例（实测） | 源码依据 |
| --- | --- | --- |
| **接触是否被焊接** | 面接触并集 `12 tri / 1 comp / vol 16`；边接触并集 `24 tri / 2 comp / vol 16`。体积相同，面积分别为 40 与 48。 | `boolean_result.cpp:142-148,284-302` |
| **共面接缝是否保留** | 面接触并集：无溯源时 `12 tri`；`AsOriginal` 时 `20 tri`（`test/boolean_test.cpp:681-689` 的 `{{12,20}}` 就是后者）。体积/面积相同（16 / 40）。 | `edge_op.cpp:239-244`；`manifold.cpp:514-528` |
| **特征是否曾经存在** | 亚 ε 薄腔消解后与实心体**逐字段相同**（`tri/area/vol/comps` 全同）。 | `edge_op.cpp:79-131` |
| **拓扑类型** | 同体积、同面积，但 1 comp vs 2 comp；genus 0 vs −1。 | `manifold.cpp:423-426` |
| **是否站在 ε 边缘** | `t=1.0·ε` 与 `t=1.2·ε` 数值相近，但一个消解、一个保留。 | `edge_op.cpp:89-96` |
| **溯源正确性** | 两个结果 `vol/tri` 相同，但一个保留了 Q 的 `faceID`，另一个把 Q 整体丢弃。 | `boolean_result.cpp:505-541` |
| **体积符号/定向** | `abs(vol)` 相同但方向相反（Subtract 时 Q 翻面），面数相同。 | `boolean_result.cpp:175-182,787` |
| **查询 ε 的可靠性** | 同一坐标在分量 ε=5e-13 下 `w=0`，在整体 ε=1e-6 下 `w=2`。 | `boolean3.cpp:747,757,777` |
| **整体是否被静默清空** | 1e12 平移下 1 单位立方体 `a+b`、`a+a` 均空且 `NoError`。 | `properties.cpp:404-407`；`boolean_result.cpp:831-833` |
| **参数是否真的生效** | `RemoveDegenerates(firstNewVert)` 的参数未实现，重复布尔会重扫旧顶点。 | `edge_op.cpp:72-79`（参数未被引用） |

---

## 10. 判断 → 源码依据 → 保证边界（逐条汇总）

| # | 判断 | 源码依据 | 保证边界（结论在什么条件下才成立） |
| --- | --- | --- | --- |
| 1 | 接触「相等」由 `min(epsP,epsQ)` 定义 | `boolean3.cpp:35-37`；`boolean2_predicates.cpp:50-55` | 仅当两个操作数的 `epsilon_` 都与其几何精度匹配时；精度错配（§7.2-L3）时该定义失去意义 |
| 2 | tie 用固定方向符号扰动打破，方向随 op 变化 | `boolean3.cpp:99-139,167-170,653-661` | `tol > 0`；`Displacements` 在 3 个以上不同面法向处退化为顶点法向（`boolean3.cpp:79-87`） |
| 3 | 面重合在 Add 下内部面消失，但接缝是否保留取决于溯源 | `boolean_result.cpp:142-148`；`test/boolean_test.cpp:681-689`（`AsOriginal` → 20 tri） | 内部面消失是几何结论（体积/面积恒定）；三角形数 12 vs 20 由 `Continuous`（`edge_op.cpp:239-244`）与 `faceID` 决定，**不是**几何结论 |
| 4 | 边/点接触在三种操作下都不焊接 | `boolean_result.cpp:284-302`；`edge_op.cpp:1119-1222`（pinch 处理） | 3D 中成立；2D 中「不焊接」表现为「输出 2 条轮廓」，而不是保留一条悬空边（`docs/Boolean2.md:135-139`） |
| 5 | Intersect 对完全重合返回 P（而非空），对零测度面接触返回空 | `boolean3.cpp:659-661`；实测 `P^P = P`、面接触 `P^Q = 空` | 仅当 `P`、`Q` 的 `inclusion` 在扰动后仍判为重叠；源码注释称「Intersect 本应收缩 Q」是刻意的 Split 性能取舍 |
| 6 | 拓扑流形 ≠ 几何接触不变 | `properties.cpp:102-128,138-140`；`boolean_result.cpp:284-302`；`edge_op.cpp:58-78` | 流形性只在组合层面成立；`IsManifold` 不看坐标；自交检查默认关闭 |
| 7 | 空结果与薄特征消解都是 `NoError` | `impl.h:145`；`boolean_result.cpp:831-833`；`properties.cpp:404-407` | 唯一可区分的错误是 `ResultTooLarge`/`Cancelled`/输入错误（`boolean_result.cpp:738-747,781-785`） |
| 8 | ε 随绝对坐标线性增长，平移不协变 | `properties.cpp:408`；`common.h:323-329`；`boolean2.cpp:272-283` | 仅当 `|T| ≲ |x|`；超出必须归一化 |
| 9 | 非均匀缩放改变特征/ε 比 | 同上（`Scale()` 取三轴最大） | 实测 `t/ε` 可从 20 降到 0.2 |
| 10 | 链式布尔会追溯性消解历史特征 | `boolean_result.cpp:933,955`；`edge_op.cpp:79-131` | 由结果 bbox 决定；且 `firstNewVert` 未实现使防护失效 |
| 11 | 相对精度不可低于 1e-12 | `utils.h:40`；`manifold.cpp:57,69-71` | 只能放松；收紧必须靠坐标归一化 |
| 12 | 2D 未做局部原点重定位（与文档不符） | `boolean2.cpp:174-177` vs `docs/Boolean2.md:149-152` | 大坐标下 2D 误差可超出舍入（实测 T=1e9 共享边并集面积 2.5 ≠ 2） |
| 13 | `MinGap()==0` 不能证明接触 | `properties.cpp:477-495` | 只可用于证明分离（`>0`） |
| 14 | 完整布尔路径上 `OriginalID` 恒为 -1；空操作数早退路径的行为不稳定 | `boolean_result.cpp:505-541,749-759`；`manifold.cpp:449-450` | 溯源审计必须用 `faceID`/`triRef`；`OriginalID` 的保留与否随空体来源与求值时机变化（§5.2） |
| 15 | debug 断言在生产构建中不存在 | `optional_assert.h:39-68` | 「manifold 输出」在 Release 下是**声明**，不是被检查的不变量；只有调用 `Is2Manifold()`/`NumDegenerateTris()` 才真正验证 |

---

## 11. 复现与验证建议

仓库内已有的、与本主题直接相关的测试（可作为验收判据的现成基线）：

| 测试 | 文件:行 | 固定的语义 |
| --- | --- | --- |
| `Boolean.FaceUnion` | `test/boolean_test.cpp:681` | 面接触 Add：1 comp、area 10、vol 2；**`AsOriginal` 下 12 vert / 20 tri** |
| `Boolean.EdgeUnion` / `EdgeUnion2` / `CornerUnion` | `:692` / `:698` / `:706` | 边或点接触 Add → 2 comps（每分量 8 vert / 12 tri） |
| `Boolean.SelfSubtract` | `:474` | 完全重合 Subtract → 空 |
| `Boolean.Coplanar` / `MultiCoplanar` | `:645` / `:660` | 共面切削的 genus 与体积 |
| `Boolean.PartialCoincidentFaceUnion` | `:1389` | 内部共享面必须删除（面积减去 2×overlap） |
| `Boolean.Precision` / `Precision2` | `:923` / `:941` | ε 作为「无交」阈值（`Precision2` 明确写 `epsilon/2` 判无交、`4·epsilon` 判有交，`test/boolean_test.cpp:941-953`） |
| `Boolean.Precision3` | `:955` | 单精度输入破坏共面性；`SetRelativePrecision(1e-7)` 恢复 |
| `Boolean.PointWithinToleranceOfCoincidentEdge` | `:1415` | ε/4 内的点判为内部，4ε 外判为外部 |
| `BooleanComplex.Close` | `test/boolean_complex_test.cpp:210` | 亚 ε 平移的 10 次自交仍保持球体积 |
| `CrossSectionComplex.CleanInputsKeepEpsSeparatedVerts` | `test/cross_section_complex_test.cpp:373` | 干净输入不得产生 ε 内的输出顶点 |
| `CrossSectionComplex.ValidatorRejects*` | `:731,742,753` | 用独立 oracle 检查保留图的边平衡 |

**建议补充的回归（对应本报告 §9 中当前实现会失败的项）：**

1. **C4 链式无侵蚀**：`fine(带 1e-9 缺口) + 远处方块` 后，再做一次结果 bbox 横跨两者的布尔，近处分量的 `NumTri/Area` 必须不变。当前失败（§7.2-L4）。
2. **B1 + R3 溯源**：`t` 从 `4ε` 扫到 `0.5ε`，断言 `faceID(切削体)` 的出现/消失只发生在 `[0.9ε, 1.2ε]` 内且只发生一次。
3. **C2 平移协变受限 + A1 状态完整**：`T ∈ {0, 1e3, 1e6}` 上同一间隙配置，断言 `|T| ≤ |x|` 时结果一致；`T = 1e9` 下 `gap=1e-4` 必须保持「分离」而不是因 ε 变大被焊接；`T = 1e12` 下 `a+a` 必须报错而不是返回 `NoError` 空体。
4. **E1/E2 消解标注**：为 `IsEmpty()` 增加「是否因 ε 消解而空」的可查询位（当前 API 缺口）。
5. **`firstNewVert` 语义**：要么实现它，要么从接口与注释中删除（当前是文档谎言，`edge_op.cpp:72-79`）。
6. **`docs/Boolean2.md` 与 `boolean2.cpp:174-177`** 的一致性：二者目前矛盾，必须择一修正。

---

## 12. 一句话回答

> Manifold 的拓扑流形保证是**组合**层面的，它保证输出的半边可以配对、法向一致、可三角化；但它**不**编码、也**不**检查任何几何接触关系。接触是否被焊接、薄特征是否被消解，由 `TieTol = min(epsP, epsQ)` 的符号扰动、`ε = relPrec · max|coord|` 的绝对坐标尺度、以及结果 bbox 驱动的 `RemoveDegenerates` 三者共同决定，而这三者对平移、缩放、链式顺序都敏感。因此验收必须是「ε 带外的组合不变量 + ε 带内的单调过渡 + 溯源集合 + 归一化后的协变」这一组关系判据；体积相等与面数相等既不能证明接触正确，也不能发现实体丢失。
