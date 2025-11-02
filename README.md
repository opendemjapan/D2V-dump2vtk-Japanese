# dump2vtk

> Single-file **GUI + CLI** tool (PySide6) for LIGGGHTS/LAMMPS dump workflows.  
> 粒子ダンプ/局所ダンプ（force chain）を **VTK/VTU** に高速変換。ヘッダー置換も同梱。

<p align="center">
<b>Tabs</b>: 1) Dump → VTK (Particles) / 2) Header replacement / 3) Dump → VTK (Force chain + Louvain)
</p>

---

## Highlights / 主要機能

- **3ツールを1つに統合（タブ分け）**
  1. Dump → VTK (粒子) — 旧 lpp.py 相当（**並列化**・チャンク・スキップ生成対応）  
  2. ヘッダー変換 — 旧 dump_rename.py 相当（`ITEM: ENTRIES ...` 置換）  
  3. Dump → VTK (force chain) — 旧 dump2force.py 相当（VTU/VTK出力＋**Louvain**＋**コミュニティ集約**）
- **巨大ファイル対応**: force chain で **ストリーミング読込（増分読込）**を実装。低メモリで処理。
- 粒子/force chain ともに **並列・チャンク処理**に対応（同じ UX）。
- **連番自動収集の改良**: `forcechain-16000-renamed.dmp` など「数字＋固定サフィックス」にも対応。
- バグ修正・堅牢化: 末尾でない数字の連番展開、NaN/Inf 処理、Windows mp 安定化 ほか。

参考: 本ツールは **Pizza.py/LPP** の考え方を継承・拡張しています。Pizza.py は **GPL-2.0**、LPP は Pizza.py の派生です。詳細はライセンス欄を参照してください。

---

## Install / 動作環境

- Python **3.8+**（推奨 3.10 以降）
- 必須: `numpy`  
- GUI を使う場合: `PySide6`  
- 任意（VTK ライブラリを使って書き出す場合）: `vtk`  
- 任意（force chain のコミュニティ検出）: `networkx`, `python-louvain`

```bash
# minimal
pip install numpy

# GUI
pip install PySide6

# optional backends / features
pip install vtk networkx python-louvain
```

---

## Quick start / 使い方

### GUI
```bash
python dump2vtk.py
```
- ドラッグ&ドロップでファイル追加、詳細設定は「Advanced」から。

### CLI
3つのサブコマンドを用意しています。

#### 1) 粒子ダンプ → VTK（VTK legacy または python-vtk）
```bash
python dump2vtk.py lpp DUMP_*.dmp \
  --format ascii --backend legacy \
  --cpunum 8 --chunksize 8 --no-overwrite
```
- `--backend vtk` を選ぶと python-vtk バックエンドを使用します（任意）。
- 出力は点群（POLYDATA）＋バウンディングボックス。

#### 2) ヘッダー変換（`ITEM: ENTRIES ...` 置換）
```bash
python dump2vtk.py rename \
  -H "ITEM: ENTRIES x1 y1 z1 x2 y2 z2 id1 id2 periodic fx fy fz" \
  forcechain-*.dmp --inplace
```

#### 3) force chain → VTK/VTU（**Louvain** + コミュニティ集約）
```bash
# VTU (XML) + binary + コミュニティ PointData も出力
python dump2vtk.py force forcechain-*-renamed.dmp \
  --vtk-format vtu --encoding binary \
  --resolution 1.0 --seed 42 --write-pointdata \
  --cpunum 8 --chunksize 8 --no-overwrite
```
- **ストリーミング読込**なので巨大ファイルでも逐次処理します。
- 1つのファイルを選ぶだけで同プレフィクスの**連番**を自動収集します。

> Louvain を使うには `networkx` と `python-louvain` が必要です。

---

## 入出力の概要

- 粒子ダンプ: `x/y/z`（または `xs/ys/zs`）を自動判定。スカラー/ベクトル列をヒューリスティックで抽出し、VTK の **SCALARS/VECTORS** として出力します。  
- force chain（pair/local dump）: 必須 12 列（`x1 y1 z1 x2 y2 z2 id1 id2 periodic fx fy fz`）＋任意列を解釈し、**線分ポリライン**として VTK/VTU に出力します。`--keep-periodic` で周期接触も保持可能。  
- Louvain でコミュニティを検出し、`community, intra_comm` を **CELL_DATA** に付与。オプションでノード側の `node_community, node_degree, node_force_sum` を **POINT_DATA** に付与。さらに各セル属性について**コミュニティ集約（mean/sum）**も計算して書き出します。

---

## Examples / 例

- **compaction_LIGGGHTS** のケースを参考にした force chain 例（ヘッダー置換 → 変換）:
```bash
# 1) ENTRIES を合わせる（必要に応じて）
python dump2vtk.py rename \
  -H "ITEM: ENTRIES x1 y1 z1 x2 y2 z2 id1 id2 periodic fx fy fz" \
  forcechain-*.dmp --inplace

# 2) 変換（VTU / binary / コミュニティ出力）
python dump2vtk.py force forcechain-*.dmp \
  --vtk-format vtu --encoding binary --write-pointdata
```

---

## Build (Windows, one-file EXE via Nuitka)

以下を PowerShell から実行（Nuitka + PySide6 プラグイン）。

```powershell
$J = [Environment]::ProcessorCount
python -m nuitka --assume-yes-for-downloads --jobs=$J --msvc=latest `
  --mode=onefile --enable-plugin=pyside6 `
  --windows-icon-from-ico=D2V.ico `
  --include-data-files="*.png=assets/" `
  --include-data-files="D2V.ico=assets/" `
  --output-dir="." --output-filename="dump2vtk.exe" --remove-output dump2vtk.py `
  --windows-console-mode=attach
```

## Windows インストーラー（Inno Setup）

```powershell
# 事前に inno setup compiler をインストール
PowerShell -NoProfile -ExecutionPolicy Bypass -File .\1_make_selfsigned_and_export_ascii.ps1 -PfxPass "powder_hunatai"
PowerShell -NoProfile -ExecutionPolicy Bypass -File .\2_build_installer_ascii_v3.ps1 -PfxPass "powder_hunatai"
```

---

## Dependencies / 依存関係

- **必須**: numpy  
- **GUI**: PySide6  
- **任意**: vtk（python-vtk バックエンド）, networkx + python-louvain（コミュニティ検出）

---

## Attribution / 謝辞

- 本ツールは **Pizza.py** の考え方（ダンプ処理・VTK 出力）および **LPP**（並列化・チャンク処理）を強く参考にしています。  
- 例は **[compaction_LIGGGHTS]** を参照しています（MIT License）。
- LIGGGHTS®/CFDEM® は各社の登録商標です。

詳しいライセンスは [`LICENCE`](./LICENCE) をご覧ください。

---

## License / ライセンス

本リポジトリは、Pizza.py 派生コードを含むため **GPL-2.0** に基づいて配布されます。  
新規に作成したコード部分については **MIT License** も併記します（ただし本リポジトリとして配布する場合は GPL-2.0 の条件が優先されます）。

SPDX: `GPL-2.0-only AND MIT`

See [`LICENCE`](./LICENCE) for the full terms.

---

## Citation (optional)

If this tool helps your work, please cite this repository as:

```
@software{dump2vtk,
  title = {dump2vtk: GUI+CLI tool for LIGGGHTS/LAMMPS dump to VTK/VTU},
  author = {dump2vtk contributors},
  year = {2025},
  url = {https://github.com/YOUR_ORG/dump2vtk}
}
```

Replace `YOUR_ORG` with your GitHub namespace when you publish.
