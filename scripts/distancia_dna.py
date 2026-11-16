"""Distância mínima (átomos pesados) de cada resíduo da proteína ao DNA, nos modelos WT do AlphaFold Server.

Uso: python -I distancia_dna.py <lote.zip> <resultados_afserver.csv>
Acrescenta ao CSV: dist_dna_wt (Å) do resíduo mutado na estrutura WT, pLDDT WT do resíduo e a classe
  'contato DNA' (<4 Å), 'próximo' (4–8 Å) ou 'distante' (>8 Å).
"""
import csv
import json
import math
import re
import sys
import zipfile
from statistics import mean

INI = {"PAX9": 6, "MSX1": 166}
WT_DIR = {"PAX9": "pax9_paired_domain_dna_wt", "MSX1": "msx1_wt"}


def ler_modelo(z, pasta):
    melhores = []
    for k in range(5):
        n = f"{pasta}/fold_{pasta}_summary_confidences_{k}.json"
        melhores.append((json.loads(z.read(n))["ranking_score"], k))
    k = max(melhores)[1]
    cif = z.read(f"{pasta}/fold_{pasta}_model_{k}.cif").decode()
    cols, prot, dna = [], {}, []
    for ln in cif.splitlines():
        if ln.startswith("_atom_site."):
            cols.append(ln.split(".", 1)[1].strip())
        elif ln.startswith(("ATOM", "HETATM")) and cols:
            r = dict(zip(cols, ln.split()))
            if r.get("type_symbol") == "H":
                continue
            xyz = (float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]))
            if r["label_asym_id"] == "A":
                prot.setdefault(int(r["label_seq_id"]), []).append((xyz, float(r["B_iso_or_equiv"])))
            else:
                dna.append(xyz)
    return prot, dna


def main():
    zpath, res_csv = sys.argv[1:3]
    z = zipfile.ZipFile(zpath)
    estr = {g: ler_modelo(z, p) for g, p in WT_DIR.items()}
    with open(res_csv, encoding="utf-8-sig") as f:
        linhas = list(csv.DictReader(f))
    for r in linhas:
        r["dist_dna_wt"] = r["plddt_wt_sitio"] = r["classe_sitio"] = ""
        if r["variante"] == "WT":
            continue
        pos = int(re.match(r"[A-Z](\d+)[A-Z]", r["variante"]).group(1))
        prot, dna = estr[r["gene"]]
        atomos = prot[pos - INI[r["gene"]] + 1]
        d = min(math.dist(a, b) for a, _ in atomos for b in dna)
        r["dist_dna_wt"] = round(d, 1)
        r["plddt_wt_sitio"] = round(mean(p for _, p in atomos), 1)
        r["classe_sitio"] = "contato DNA" if d < 4 else ("próximo" if d < 8 else "distante")
    with open(res_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()))
        w.writeheader()
        w.writerows(linhas)
    for r in linhas:
        print(r["job"], r["dist_dna_wt"], r["classe_sitio"], r["plddt_wt_sitio"])


if __name__ == "__main__":
    main()
