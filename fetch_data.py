"""
Pull Allen Mouse Brain Atlas ISH expression data for presynaptic / dopamine
pathway genes across brain regions relevant to Parkinson's disease and
the Edwards lab's research on VMAT2, alpha-synuclein, and GBA.

Genes:
  Core PD/presynaptic:  Snca, Slc18a2 (VMAT2), Th, Gba, Slc6a3 (DAT)
  Corelease:            Slc17a6 (VGLUT2), Slc17a7 (VGLUT1)
  Synuclein family:     Sncb, Sncg
  Dense core / peptide: Pcsk1, Cpe

Structures:
  Substantia nigra (SN) 374, VTA 749, Caudoputamen (CP) 672,
  Nucleus accumbens (ACB) 56, Isocortex 315, Hippocampus 1089,
  Locus coeruleus (LC) 372, Dorsal raphe (DR) 872
"""
import requests, csv, time
from pathlib import Path

OUT  = Path(__file__).parent
BASE = "http://api.brain-map.org/api/v2/data/query.json"

GENES = {
    "core": ["Snca", "Slc18a2", "Th", "Gba", "Slc6a3"],
    "corelease": ["Slc17a6", "Slc17a7"],
    "synuclein": ["Sncb", "Sncg"],
    "dcv": ["Pcsk1", "Cpe"],
}
ALL_GENES = [g for grp in GENES.values() for g in grp]

STRUCTURES = {
    "Subst_nigra":  374,
    "VTA":          749,
    "Caudoputamen": 672,
    "Nuc_accumbens": 56,
    "Isocortex":    315,
    "Hippocampus":  1089,
    "Locus_coer":   147,
    "Dorsal_raphe": 872,
}

def get_dataset(gene):
    url = (f"{BASE}?criteria=model::SectionDataSet,"
           f"rma::criteria,genes[acronym$eq'{gene}'],"
           f"[failed$eqfalse],"
           f"products[abbreviation$eqMouse]"
           f"&num_rows=5")
    r = requests.get(url, timeout=20).json()
    if not r.get("msg"):
        return None
    for row in r["msg"]:
        if row.get("plane_of_section_id") == 2:
            return row["id"]
    return r["msg"][0]["id"] if r["msg"] else None

def get_expression(dataset_id, struct_ids):
    id_str = ",".join(str(i) for i in struct_ids)
    url = (f"{BASE}?criteria=model::StructureUnionize,"
           f"rma::criteria,[section_data_set_id$eq{dataset_id}],"
           f"[structure_id$in{id_str}]&num_rows=50")
    r = requests.get(url, timeout=20).json()
    out = {}
    for row in r.get("msg", []):
        out[row["structure_id"]] = row.get("expression_energy", 0.0) or 0.0
    return out

struct_ids  = list(STRUCTURES.values())
id_to_name  = {v: k for k, v in STRUCTURES.items()}
gene_to_grp = {g: grp for grp, gs in GENES.items() for g in gs}

rows = []
for gene in ALL_GENES:
    print(f"  {gene}...", end=" ", flush=True)
    ds = get_dataset(gene)
    if ds is None:
        print("no dataset"); continue
    expr = get_expression(ds, struct_ids)
    row = {"gene": gene, "group": gene_to_grp[gene], "dataset_id": ds}
    for sid, name in id_to_name.items():
        row[name] = round(expr.get(sid, 0.0), 4)
    rows.append(row)
    print(f"ds={ds} ✓")
    time.sleep(0.3)

fieldnames = ["gene","group","dataset_id"] + list(STRUCTURES.keys())
with open(OUT / "expression.tsv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"\nSaved expression.tsv ({len(rows)} genes)")
for r in rows:
    vals = {k: r[k] for k in STRUCTURES}
    top  = max(vals, key=vals.get)
    print(f"  {r['gene']:12s}  peak={top} ({vals[top]:.3f})")
