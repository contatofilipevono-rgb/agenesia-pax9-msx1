"""Melhoria 5: robustez da regra (contato OU ΔΔG>=corte).

(a) Validação cruzada leave-one-out do corte de ΔΔG: em cada rodada, escolhe o corte (0 a 3, passo 0,25) que
    maximiza Youden (sens + espec - 1) nos demais e testa na variante deixada de fora. Reporta sens/espec fora da amostra.
(b) Sensibilidade ao corte de distância (3,5 / 4 / 5 Å) e de ΔΔG (0,5 / 1 / 1,5).
(c) Distância ao DNA nas 5 amostras do modelo WT (variação do contato).
Usa grupo_controle.csv (patogênicas x controle gnomAD AC>=2), só PAX9 (MSX1 tem n pequeno).
Uso: python -I -X utf8 robustez.py <grupo_controle.csv> <lote.zip>
"""
import csv
import json
import math
import sys
import zipfile


def regra(r, cd=4.0, cg=1.0):
    d = float(r["dist_dna"])
    g = float(r["ddg_thermompnn"]) if r["ddg_thermompnn"] != "" else -9
    return d < cd or g >= cg


def main():
    rows = [r for r in csv.DictReader(open(sys.argv[1], encoding="utf-8-sig")) if r["gene"] == "PAX9"]
    P = [r for r in rows if r["grupo"] == "patogenica"]
    C = [r for r in rows if r["grupo"] == "controle_gnomad" and int(r["gnomad_ac"]) >= 2]
    todos = [(r, 1) for r in P] + [(r, 0) for r in C]
    cortes = [i * 0.25 for i in range(13)]
    tp = fn = tn = fp = 0
    escolhidos = []
    for i, (r, y) in enumerate(todos):
        treino = [t for j, t in enumerate(todos) if j != i]
        def youden(c):
            pp = [regra(t, cg=c) for t, yy in treino if yy == 1]
            cc = [not regra(t, cg=c) for t, yy in treino if yy == 0]
            return sum(pp) / len(pp) + sum(cc) / len(cc) - 1
        c = max(cortes, key=youden)
        escolhidos.append(c)
        pred = regra(r, cg=c)
        if y == 1:
            tp += pred
            fn += not pred
        else:
            tn += not pred
            fp += pred
    sens, esp = tp / (tp + fn), tn / (tn + fp)
    print(f"(a) LOO: cortes escolhidos {min(escolhidos)}–{max(escolhidos)} (moda {max(set(escolhidos), key=escolhidos.count)});"
          f" fora da amostra: sensibilidade {sens:.2f}, especificidade {esp:.2f}, LR+ {sens/(1-esp):.1f}, LR- {(1-sens)/esp:.2f}")
    print("(b) sensibilidade aos cortes (sens / espec):")
    for cd in (3.5, 4.0, 5.0):
        lin = []
        for cg in (0.5, 1.0, 1.5):
            s = sum(regra(r, cd, cg) for r in P) / len(P)
            e = sum(not regra(r, cd, cg) for r in C) / len(C)
            lin.append(f"ΔΔG≥{cg}: {s:.2f}/{e:.2f}")
        print(f"    dist<{cd} Å  " + "  ".join(lin))
    # (c) distâncias nas 5 amostras
    z = zipfile.ZipFile(sys.argv[2])
    p = "pax9_paired_domain_dna_wt"
    contato_por_amostra = []
    for k in range(5):
        cif = z.read(f"{p}/fold_{p}_model_{k}.cif").decode()
        cols, prot, dna = [], {}, []
        for ln in cif.splitlines():
            if ln.startswith("_atom_site."):
                cols.append(ln.split(".", 1)[1].strip())
            elif ln.startswith(("ATOM", "HETATM")) and cols:
                rr = dict(zip(cols, ln.split()))
                if rr.get("type_symbol") == "H":
                    continue
                xyz = (float(rr["Cartn_x"]), float(rr["Cartn_y"]), float(rr["Cartn_z"]))
                if rr["label_asym_id"] == "A":
                    prot.setdefault(int(rr["label_seq_id"]), []).append(xyz)
                else:
                    dna.append(xyz)
        cont = {i + 5 for i, at in prot.items() if min(math.dist(a, b) for a in at for b in dna) < 4}
        contato_por_amostra.append(cont)
    sitios_p = {int(r["variante"][1:-1]) for r in P}
    estaveis = set.intersection(*[c & sitios_p for c in contato_por_amostra])
    variaveis = set.union(*[c & sitios_p for c in contato_por_amostra]) - estaveis
    print(f"(c) sítios patogênicos em contato em TODAS as 5 amostras: {sorted(estaveis)}; variam entre amostras: {sorted(variaveis)}")


if __name__ == "__main__":
    main()
