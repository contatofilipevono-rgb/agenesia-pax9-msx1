"""Melhoria 3: quantifica o efeito de cada mutação na interface proteína–DNA usando os modelos do AlphaFold Server.

Para cada job (melhor amostra por ranking_score) conta contatos proteína–DNA (pares de átomos pesados < 4 Å) e
resíduos da proteína na interface; compara mutante x WT do mesmo gene. Também usa o ipTM proteína–DNA e o PAE mínimo.
Uso: python -I -X utf8 interface_dna.py <lote1.zip> <lote2.zip> <resultados_afserver.csv>
"""
import csv
import json
import math
import sys
import zipfile


def modelo(z, pasta):
    k = max((json.loads(z.read(f"{pasta}/fold_{pasta}_summary_confidences_{i}.json"))["ranking_score"], i)
            for i in range(5))[1]
    cif = z.read(f"{pasta}/fold_{pasta}_model_{k}.cif").decode()
    cols, prot, dna = [], [], []
    for ln in cif.splitlines():
        if ln.startswith("_atom_site."):
            cols.append(ln.split(".", 1)[1].strip())
        elif ln.startswith(("ATOM", "HETATM")) and cols:
            r = dict(zip(cols, ln.split()))
            if r.get("type_symbol") == "H":
                continue
            xyz = (float(r["Cartn_x"]), float(r["Cartn_y"]), float(r["Cartn_z"]))
            (prot if r["label_asym_id"] == "A" else dna).append((int(r["label_seq_id"]), xyz))
    return prot, dna


def interface(prot, dna, corte=4.0):
    pares, residuos = 0, set()
    # grade simples para acelerar
    grade = {}
    for _, (x, y, zz) in dna:
        grade.setdefault((int(x // corte), int(y // corte), int(zz // corte)), []).append((x, y, zz))
    for res, (x, y, zz) in prot:
        cx, cy, cz = int(x // corte), int(y // corte), int(zz // corte)
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    for b in grade.get((cx + dx, cy + dy, cz + dz), []):
                        if math.dist((x, y, zz), b) < corte:
                            pares += 1
                            residuos.add(res)
    return pares, len(residuos)


def main():
    z1, z2 = zipfile.ZipFile(sys.argv[1]), zipfile.ZipFile(sys.argv[2])
    res_csv = sys.argv[3]
    pastas = {}
    for z in (z1, z2):
        for n in z.namelist():
            if "/" in n:
                pastas[n.split("/")[0]] = z
    rows = list(csv.DictReader(open(res_csv, encoding="utf-8-sig")))
    nome = {"PAX9_WT": "pax9_paired_domain_dna_wt"}
    base = {}
    for r in rows:
        p = nome.get(r["job"], r["job"].lower())
        if p not in pastas:
            continue
        prot, dna = modelo(pastas[p], p)
        pares, nres = interface(prot, dna)
        r["contatos_dna"] = pares
        r["residuos_interface"] = nres
        if r["variante"] == "WT":
            base[r["gene"]] = (pares, nres, float(r["iptm_prot_dna"]))
    for r in rows:
        if r.get("contatos_dna") in (None, "") or r["gene"] not in base:
            continue
        b = base[r["gene"]]
        r["delta_contatos_pct"] = round(100 * (int(r["contatos_dna"]) - b[0]) / b[0], 1)
        r["delta_iptm_dna"] = round(float(r["iptm_prot_dna"]) - b[2], 3)
    with open(res_csv, "w", newline="", encoding="utf-8-sig") as f:
        campos = list(rows[0].keys())
        for c in ("contatos_dna", "residuos_interface", "delta_contatos_pct", "delta_iptm_dna"):
            if c not in campos:
                campos.append(c)
        w = csv.DictWriter(f, fieldnames=campos)
        w.writeheader()
        w.writerows(rows)
    for g in ("PAX9", "MSX1"):
        print(f"{g} WT: {base[g][0]} pares de contato, {base[g][1]} resíduos na interface, ipTM DNA {base[g][2]}")
        sub = [r for r in rows if r["gene"] == g and r["variante"] != "WT" and r.get("delta_contatos_pct") not in (None, "")]
        for r in sorted(sub, key=lambda r: float(r["delta_contatos_pct"])):
            print(f"  {r['variante']:6s} {r['classe_sitio']:12s} Δcontatos={float(r['delta_contatos_pct']):6.1f}% "
                  f"ΔipTM={float(r['delta_iptm_dna']):+.3f} ΔΔG={r.get('ddg_thermompnn', '')}")


if __name__ == "__main__":
    main()
