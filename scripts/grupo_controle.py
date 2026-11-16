"""Grupo controle: variantes missense do gnomAD v4 dentro dos domínios de ligação ao DNA,
que NÃO são patogênicas (nem P no ClinVar nem P-lit). Compara com as patogênicas em:
distância ao DNA (estrutura WT), ΔΔG (ThermoMPNN) e AlphaMissense.

Uso: python -I -X utf8 grupo_controle.py <variantes_PAX9_MSX1.csv> <lote.zip> <pasta_estabilidade> <saida.csv>
Controle "todos": observadas no gnomAD (AC>=1). Controle "recorrentes": AC>=2 (menos chance de serem patogênicas raras).
"""
import csv
import math
import os
import runpy
import sys
import zipfile
from statistics import median

INI = {"PAX9": 6, "MSX1": 166}
FIM = {"PAX9": 130, "MSX1": 231}  # domínios (UniProt): paired 4-130, homeobox 172-231 (limitados à construção)
DOM_INI = {"PAX9": 6, "MSX1": 172}


def mw(a, b):
    """AUC (a > b) e p bicaudal aproximado (Mann-Whitney, aproximação normal)."""
    U = sum((x > y) + 0.5 * (x == y) for x in a for y in b)
    n1, n2 = len(a), len(b)
    z = (U - n1 * n2 / 2) / math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    return U / (n1 * n2), math.erfc(abs(z) / math.sqrt(2))


def fisher_2x2(a, b, c, d):
    """p bicaudal do teste exato de Fisher para [[a,b],[c,d]]."""
    from math import comb
    n = a + b + c + d
    r1, c1 = a + b, a + c
    def p(x):
        return comb(c1, x) * comb(n - c1, r1 - x) / comb(n, r1)
    p0 = p(a)
    lo, hi = max(0, r1 + c1 - n), min(r1, c1)
    return sum(p(x) for x in range(lo, hi + 1) if p(x) <= p0 + 1e-12)


def main():
    var_csv, zpath, pasta_est, saida = sys.argv[1:5]
    dist_mod = runpy.run_path(os.path.join(os.path.dirname(__file__), "distancia_dna.py"))
    z = zipfile.ZipFile(zpath)
    dist = {}
    for g, p in dist_mod["WT_DIR"].items():
        prot, dna = dist_mod["ler_modelo"](z, p)
        for i, at in prot.items():
            dist[(g, i + INI[g] - 1)] = min(math.dist(a, b) for a, _ in at for b in dna)
    ddg = {}
    for g in INI:
        path = os.path.join(pasta_est, f"ThermoMPNN_inference_{g}_WT_cadeiaA.csv")
        for r in csv.DictReader(open(path, encoding="utf-8")):
            if r["wildtype"] != r["mutation"]:
                ddg[(g, int(r["position"]) + INI[g], r["mutation"])] = float(r["ddG_pred"])  # 0-based -> numeração da proteína
    linhas = []
    for r in csv.DictReader(open(var_csv, encoding="utf-8-sig")):
        g, v = r["gene"], r["variante"]
        pos = int(r["posicao"])
        if not (DOM_INI[g] <= pos <= FIM[g]) or v.startswith("M1"):
            continue
        ac = int(r["gnomad_ac"] or 0)
        if r["rotulo"] in ("P", "P-lit"):
            grupo = "patogenica"
        elif ac >= 1:
            grupo = "controle_gnomad"
        else:
            continue
        d = dist.get((g, pos))
        dd = ddg.get((g, pos, v[-1]))
        linhas.append({"gene": g, "variante": v, "grupo": grupo, "gnomad_ac": ac, "gnomad_af": r["gnomad_af"],
                       "clinvar": r["clinvar"], "dist_dna": round(d, 1) if d is not None else "",
                       "contato_dna": "sim" if d is not None and d < 4 else "nao",
                       "ddg_thermompnn": round(dd, 2) if dd is not None else "",
                       "alphamissense": r["alphamissense"]})
    with open(saida, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()))
        w.writeheader()
        w.writerows(linhas)

    for g in INI:
        for nome, filtro in (("todos (AC>=1)", lambda r: True), ("recorrentes (AC>=2)", lambda r: r["gnomad_ac"] >= 2)):
            P = [r for r in linhas if r["gene"] == g and r["grupo"] == "patogenica"]
            C = [r for r in linhas if r["gene"] == g and r["grupo"] == "controle_gnomad" and filtro(r)]
            if len(C) < 3:
                print(f"{g} controle {nome}: só {len(C)} variantes — insuficiente")
                continue
            pc = sum(r["contato_dna"] == "sim" for r in P)
            cc = sum(r["contato_dna"] == "sim" for r in C)
            p_f = fisher_2x2(pc, len(P) - pc, cc, len(C) - cc)
            dP = [r["ddg_thermompnn"] for r in P if r["ddg_thermompnn"] != ""]
            dC = [r["ddg_thermompnn"] for r in C if r["ddg_thermompnn"] != ""]
            aP = [float(r["alphamissense"]) for r in P]
            aC = [float(r["alphamissense"]) for r in C]
            # mecanismo combinado: contato OU ΔΔG>=1
            mP = sum(r["contato_dna"] == "sim" or (r["ddg_thermompnn"] != "" and r["ddg_thermompnn"] >= 1) for r in P)
            mC = sum(r["contato_dna"] == "sim" or (r["ddg_thermompnn"] != "" and r["ddg_thermompnn"] >= 1) for r in C)
            p_m = fisher_2x2(mP, len(P) - mP, mC, len(C) - mC)
            auc_d, p_d = mw(dP, dC)
            auc_a, p_a = mw(aP, aC)
            print(f"\n{g} — patogênicas (n={len(P)}) x controle gnomAD {nome} (n={len(C)})")
            print(f"  contato com DNA: {pc}/{len(P)} ({pc/len(P):.0%}) x {cc}/{len(C)} ({cc/len(C):.0%})  Fisher p={p_f:.3g}")
            print(f"  ΔΔG mediano: {median(dP):.2f} x {median(dC):.2f}  AUC={auc_d:.2f} p={p_d:.3g}")
            print(f"  contato OU ΔΔG>=1: {mP}/{len(P)} ({mP/len(P):.0%}) x {mC}/{len(C)} ({mC/len(C):.0%})  Fisher p={p_m:.3g}")
            print(f"  AlphaMissense mediano: {median(aP):.3f} x {median(aC):.3f}  AUC={auc_a:.2f} p={p_a:.3g}")


if __name__ == "__main__":
    main()
