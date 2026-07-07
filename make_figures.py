"""
Two SVG figures for the Edwards lab demo:
1. pd_heatmap.svg   — all 11 genes × 8 regions, grouped by function, SN+VTA highlighted
2. corelease.svg    — focused bar chart: SN vs VTA for VMAT2, VGLUT2, TH, GBA, DAT, SNCA
"""
import csv, math
from pathlib import Path

OUT = Path(__file__).parent

STRUCTURES = ["Subst_nigra","VTA","Caudoputamen","Nuc_accumbens",
              "Isocortex","Hippocampus","Locus_coer","Dorsal_raphe"]
STRUCT_LABELS = ["Subst.\nnigra","VTA","Caudate-\nputamen","Nucleus\naccumbens",
                 "Isocortex","Hippocampus","Locus\ncoer.","Dorsal\nraphe"]

GROUPS = [
    ("core PD / presynaptic",   ["Snca","Slc18a2","Th","Gba","Slc6a3"],   "#c0392b"),
    ("corelease",               ["Slc17a6","Slc17a7"],                     "#2980b9"),
    ("synuclein family",        ["Sncb","Sncg"],                           "#8e44ad"),
    ("dense-core / peptide",    ["Pcsk1","Cpe"],                           "#27ae60"),
]
GROUP_ORDER = [g for _, gs, _ in GROUPS for g in gs]
GENE_COLOR  = {g: c for _, gs, c in GROUPS for g in gs}
GENE_GROUP  = {g: lbl for lbl, gs, _ in GROUPS for g in gs}

rows_by_gene = {}
with open(OUT / "expression.tsv") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        rows_by_gene[r["gene"]] = [float(r[s]) for s in STRUCTURES]


def warm_cold(v):
    if v <= 0.5:
        t = v / 0.5
        return f"rgb({int(220+(255-220)*t)},{int(235+(255-235)*t)},255)"
    else:
        t = (v - 0.5) / 0.5
        return f"rgb(255,{int(255+(70-255)*t)},{int(255+(70-255)*t)})"


def row_norm(vals):
    mn, mx = min(vals), max(vals)
    return [(v - mn) / (mx - mn) if mx > mn else 0.0 for v in vals]


# ── Figure 1: full heatmap ─────────────────────────────────────────────────────
CELL_W = 64
CELL_H = 26
LEFT   = 88
TOP    = 100
n_genes   = len(GROUP_ORDER)
n_structs = len(STRUCTURES)
W = LEFT + n_structs * CELL_W + 70
H = TOP + n_genes * CELL_H + 100

# Highlight SN and VTA columns
col_bg = ""
for si in [0, 1]:   # SN, VTA
    x = LEFT + si * CELL_W
    col_bg += (f'<rect x="{x}" y="{TOP-48}" width="{CELL_W}" height="{n_genes*CELL_H+52}" '
               f'fill="#fff8f0" rx="0" opacity="0.7"/>')

cells = ""
group_dividers = ""
prev_grp = None
gi_abs = 0
for _, gs, _ in GROUPS:
    for gene in gs:
        if gene not in rows_by_gene:
            gi_abs += 1
            continue
        vals = rows_by_gene[gene]
        nv   = row_norm(vals)
        for si, v in enumerate(nv):
            x = LEFT + si * CELL_W
            y = TOP  + gi_abs * CELL_H
            color = warm_cold(v)
            cells += (f'<rect x="{x}" y="{y}" width="{CELL_W}" height="{CELL_H}" '
                      f'fill="{color}" stroke="white" stroke-width="1.5"/>')
            raw = vals[si]
            if raw > 1.0 or (v > 0.85 and raw > 0.05):
                fc = "white" if v > 0.88 else "#333"
                cells += (f'<text x="{x+CELL_W/2:.0f}" y="{y+CELL_H/2+4:.0f}" '
                          f'text-anchor="middle" font-size="8" fill="{fc}">{raw:.2f}</text>')
        gi_abs += 1
    # group divider line
    group_dividers += (f'<line x1="{LEFT}" y1="{TOP+gi_abs*CELL_H}" '
                       f'x2="{LEFT+n_structs*CELL_W}" y2="{TOP+gi_abs*CELL_H}" '
                       f'stroke="#ccc" stroke-width="1.5"/>')

# Row labels
row_lbls = ""
gi_abs = 0
for _, gs, grp_color in GROUPS:
    # group label on left
    mid_y = TOP + gi_abs * CELL_H + len(gs) * CELL_H / 2 + 4
    row_lbls += (f'<text x="4" y="{mid_y:.0f}" font-size="9" fill="{grp_color}" '
                 f'font-weight="600" writing-mode="horizontal-tb">'
                 f'{GENE_GROUP[gs[0]]}</text>')
    for gene in gs:
        y = TOP + gi_abs * CELL_H + CELL_H // 2 + 4
        row_lbls += (f'<text x="{LEFT-5}" y="{y}" text-anchor="end" font-size="10.5" '
                     f'fill="{GENE_COLOR[gene]}" font-weight="600" '
                     f'font-family="ui-monospace,Menlo,monospace">{gene}</text>')
        gi_abs += 1

# Column headers
col_hdrs = ""
for si, lbl in enumerate(STRUCT_LABELS):
    x = LEFT + si * CELL_W + CELL_W // 2
    parts = lbl.split("\n")
    for li, part in enumerate(parts):
        y = TOP - 10 - (len(parts) - 1 - li) * 13
        fc = "#c05000" if si in (0, 1) else "#444"
        fw = "bold"    if si in (0, 1) else "normal"
        col_hdrs += (f'<text x="{x}" y="{y}" text-anchor="middle" font-size="9.5" '
                     f'fill="{fc}" font-weight="{fw}">{part}</text>')

# Colorbar
cb_x = LEFT + n_structs * CELL_W + 10
cb_h = n_genes * CELL_H
colorbar = (f'<defs><linearGradient id="cb1" x1="0" y1="1" x2="0" y2="0">'
            f'<stop offset="0%" stop-color="{warm_cold(0)}"/>'
            f'<stop offset="50%" stop-color="{warm_cold(0.5)}"/>'
            f'<stop offset="100%" stop-color="{warm_cold(1)}"/>'
            f'</linearGradient></defs>'
            f'<rect x="{cb_x}" y="{TOP}" width="13" height="{cb_h}" '
            f'fill="url(#cb1)" stroke="#ccc" stroke-width="0.5"/>'
            f'<text x="{cb_x+6}" y="{TOP-4}" text-anchor="middle" font-size="8" fill="#555">high</text>'
            f'<text x="{cb_x+6}" y="{TOP+cb_h+10}" text-anchor="middle" font-size="8" fill="#555">low</text>')

svg1 = f"""<svg viewBox="0 0 {W} {H}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{W//2}" y="22" text-anchor="middle" font-size="13" font-weight="600" fill="#222">
    Presynaptic and PD-Relevant Gene Expression Across Brain Regions
  </text>
  <text x="{W//2}" y="38" text-anchor="middle" font-size="10" fill="#666">
    Allen Mouse Brain Atlas · ISH expression energy (row-normalized per gene)
  </text>
  <text x="{W//2}" y="53" text-anchor="middle" font-size="10" fill="#c05000">
    ▶ Substantia nigra and VTA highlighted — primary sites of PD-relevant dopamine neuron pathology
  </text>
  {col_bg}
  {colorbar}
  {col_hdrs}
  {cells}
  {group_dividers}
  {row_lbls}
</svg>"""

with open(OUT / "pd_heatmap.svg", "w") as f:
    f.write(svg1)
print("Wrote pd_heatmap.svg")


# ── Figure 2: SN vs VTA focused bar chart ────────────────────────────────────
# Show raw expression for 6 key genes in SN and VTA side by side
FOCUS_GENES  = ["Th", "Gba", "Slc18a2", "Slc6a3", "Snca", "Slc17a6"]
FOCUS_LABELS = ["Th\n(TH)","Gba\n(GBA)","Slc18a2\n(VMAT2)","Slc6a3\n(DAT)","Snca\n(α-syn)","Slc17a6\n(VGLUT2)"]
FOCUS_COLORS = ["#c0392b","#e67e22","#2980b9","#16a085","#8e44ad","#27ae60"]

SN_IDX  = STRUCTURES.index("Subst_nigra")
VTA_IDX = STRUCTURES.index("VTA")

sn_vals  = [rows_by_gene[g][SN_IDX]  for g in FOCUS_GENES]
vta_vals = [rows_by_gene[g][VTA_IDX] for g in FOCUS_GENES]

FW2 = 660
FH2 = 400
PAD_L = 55
PAD_B = 80
PAD_T = 80
AREA_W2 = FW2 - PAD_L - 30
AREA_H2 = FH2 - PAD_B - PAD_T

n_groups = len(FOCUS_GENES)
GROUP_W  = AREA_W2 / n_groups
BAR_W    = GROUP_W * 0.32
GAP      = GROUP_W * 0.04

max_v2 = max(max(sn_vals), max(vta_vals)) * 1.12
y_scale2 = AREA_H2 / max_v2

def bx(gi): return PAD_L + gi * GROUP_W + GROUP_W * 0.1
def by(v):  return PAD_T + AREA_H2 - v * y_scale2

bars2 = ""
for i, (gene, lbl, col) in enumerate(zip(FOCUS_GENES, FOCUS_LABELS, FOCUS_COLORS)):
    x_sn  = bx(i)
    x_vta = x_sn + BAR_W + GAP
    h_sn  = sn_vals[i]  * y_scale2
    h_vta = vta_vals[i] * y_scale2
    # SN bar (solid)
    bars2 += (f'<rect x="{x_sn:.1f}" y="{by(sn_vals[i]):.1f}" width="{BAR_W:.1f}" '
              f'height="{h_sn:.1f}" fill="{col}" opacity="0.9" rx="2"/>')
    # VTA bar (hatched feel — lighter)
    bars2 += (f'<rect x="{x_vta:.1f}" y="{by(vta_vals[i]):.1f}" width="{BAR_W:.1f}" '
              f'height="{h_vta:.1f}" fill="{col}" opacity="0.45" rx="2" '
              f'stroke="{col}" stroke-width="1.5"/>')
    # value labels
    if sn_vals[i] > 0.01:
        bars2 += (f'<text x="{x_sn+BAR_W/2:.1f}" y="{by(sn_vals[i])-4:.1f}" '
                  f'text-anchor="middle" font-size="9" fill="{col}">{sn_vals[i]:.2f}</text>')
    if vta_vals[i] > 0.01:
        bars2 += (f'<text x="{x_vta+BAR_W/2:.1f}" y="{by(vta_vals[i])-4:.1f}" '
                  f'text-anchor="middle" font-size="9" fill="{col}">{vta_vals[i]:.2f}</text>')
    # gene label
    parts = lbl.split("\n")
    for li, part in enumerate(parts):
        ly = PAD_T + AREA_H2 + 14 + li * 13
        bars2 += (f'<text x="{x_sn+BAR_W+GAP/2:.1f}" y="{ly}" text-anchor="middle" '
                  f'font-size="10" fill="{col}" font-weight="600" '
                  f'font-family="ui-monospace,Menlo,monospace">{part}</text>')

# Y-axis ticks
yticks2 = ""
for t in [0, 1, 2, 3, 4, 5, 6, 7]:
    if t > max_v2: break
    yy = by(t)
    yticks2 += (f'<line x1="{PAD_L-4}" y1="{yy:.1f}" x2="{PAD_L+AREA_W2}" y2="{yy:.1f}" '
                f'stroke="#eee" stroke-width="1"/>'
                f'<text x="{PAD_L-7}" y="{yy+3:.1f}" text-anchor="end" font-size="9" fill="#666">{t}</text>')

# Legend
legend = (f'<rect x="{PAD_L}" y="58" width="12" height="10" fill="#555" opacity="0.9" rx="2"/>'
          f'<text x="{PAD_L+16}" y="67" font-size="10" fill="#333">Substantia nigra</text>'
          f'<rect x="{PAD_L+130}" y="58" width="12" height="10" fill="#555" opacity="0.45" '
          f'stroke="#555" stroke-width="1.5" rx="2"/>'
          f'<text x="{PAD_L+146}" y="67" font-size="10" fill="#333">VTA</text>')

# VGLUT2 annotation arrow
vglut_i = FOCUS_GENES.index("Slc17a6")
ann_x   = bx(vglut_i) + BAR_W + GAP / 2
ann_y   = by(vta_vals[vglut_i]) - 22
bars2 += (f'<text x="{ann_x:.1f}" y="{ann_y:.1f}" text-anchor="middle" '
          f'font-size="9" fill="#27ae60" font-style="italic">'
          f'VTA &gt; SN:</text>'
          f'<text x="{ann_x:.1f}" y="{ann_y+12:.1f}" text-anchor="middle" '
          f'font-size="9" fill="#27ae60" font-style="italic">corelease</text>')

svg2 = f"""<svg viewBox="0 0 {FW2} {FH2}" xmlns="http://www.w3.org/2000/svg"
     style="font-family:-apple-system,system-ui,sans-serif;background:white;">
  <text x="{FW2//2}" y="22" text-anchor="middle" font-size="13" font-weight="600" fill="#222">
    Key PD Genes: Substantia Nigra vs VTA Expression
  </text>
  <text x="{FW2//2}" y="38" text-anchor="middle" font-size="10" fill="#666">
    Raw ISH expression energy · VGLUT2 enriched in VTA captures dopamine–glutamate corelease specificity
  </text>
  {legend}
  <line x1="{PAD_L}" y1="{PAD_T}" x2="{PAD_L}" y2="{PAD_T+AREA_H2}" stroke="#aaa" stroke-width="1.5"/>
  <line x1="{PAD_L}" y1="{PAD_T+AREA_H2}" x2="{PAD_L+AREA_W2}" y2="{PAD_T+AREA_H2}" stroke="#aaa" stroke-width="1.5"/>
  {yticks2}
  {bars2}
  <text transform="rotate(-90)" x="-{PAD_T+AREA_H2//2}" y="14" text-anchor="middle"
        font-size="10.5" fill="#444">ISH expression energy</text>
</svg>"""

with open(OUT / "corelease.svg", "w") as f:
    f.write(svg2)
print("Wrote corelease.svg")
