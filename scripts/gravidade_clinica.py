"""Melhoria 2: mecanismo previsto x gravidade clínica (nº de dentes ausentes, sem terceiros molares).

Fonte: Suplemento de Intarak et al. 2023 (tabelas já extraídas em intarak2023_supl.txt). Usamos, para cada variante,
o nº de dentes ausentes da linha da variante (paciente índice ou valor relatado). Também comparamos com variantes
nulas (frameshift/nonsense = perda de uma cópia), que servem de referência de gravidade máxima por haploinsuficiência.
Uso: python -I -X utf8 gravidade_clinica.py <intarak2023_supl.txt> <resultados_afserver.csv>
"""
import csv
import math
import re
import sys
from statistics import median

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H",
       "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
       "Tyr": "Y", "Val": "V"}


def mw(a, b):
    U = sum((x > y) + 0.5 * (x == y) for x in a for y in b)
    n1, n2 = len(a), len(b)
    z = (U - n1 * n2 / 2) / math.sqrt(n1 * n2 * (n1 + n2 + 1) / 12)
    return U / (n1 * n2), math.erfc(abs(z) / math.sqrt(2))


def main():
    txt, res_csv = sys.argv[1:3]
    corte = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0  # corte de ΔΔG (artigo: 1,5)
    missense, nulas = {}, []
    for ln in open(txt, encoding="utf-8"):
        cols = [c.strip() for c in ln.split("|")]
        if len(cols) < 6:
            continue
        aa = cols[2]
        try:
            dentes = float(cols[4])
        except ValueError:
            continue
        m = re.fullmatch(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", aa)
        if m and m.group(1) in AA3 and m.group(3) in AA3:
            v = AA3[m.group(1)] + m.group(2) + AA3[m.group(3)]
            missense.setdefault(v, dentes)
        elif re.search(r"fs\*|\*$|Ter", aa) or "fs" in aa:
            nulas.append(dentes)
    res = {r["variante"]: r for r in csv.DictReader(open(res_csv, encoding="utf-8-sig")) if r["gene"] == "PAX9"}
    grupos = {"contato DNA": [], "desestabiliza (ΔΔG≥1, sem contato)": [], "nenhum dos dois": []}
    linhas = []
    for v, d in missense.items():
        r = res.get(v)
        if not r or r.get("ddg_thermompnn") in ("", None):
            continue
        ddg = float(r["ddg_thermompnn"])
        if r["classe_sitio"] == "contato DNA":
            g = "contato DNA"
        elif ddg >= corte:
            g = "desestabiliza (ΔΔG≥1, sem contato)"
        else:
            g = "nenhum dos dois"
        grupos[g].append(d)
        linhas.append((v, g, ddg, d))
    for v, g, ddg, d in sorted(linhas, key=lambda x: -x[3]):
        print(f"  {v:6s} {g:36s} ΔΔG={ddg:5.2f} dentes={d:g}")
    print("\nDentes ausentes (mediana, n):")
    for g, vals in grupos.items():
        if vals:
            print(f"  {g:36s} {median(vals):5.1f}  (n={len(vals)})")
    print(f"  {'variantes nulas (frameshift/nonsense)':36s} {median(nulas):5.1f}  (n={len(nulas)})")
    a, b = grupos["contato DNA"], grupos["desestabiliza (ΔΔG≥1, sem contato)"]
    if a and b:
        auc, p = mw(b, a)
        print(f"\nDesestabilizantes x contato: AUC={auc:.2f}, p={p:.3f}")
    todas_mis = [d for vals in grupos.values() for d in vals]
    auc, p = mw(nulas, todas_mis)
    print(f"Nulas x missense: AUC={auc:.2f}, p={p:.3f}")
    # correlação ΔΔG x dentes
    xs = [x[2] for x in linhas]
    ys = [x[3] for x in linhas]
    def rank(x):
        o = sorted(range(len(x)), key=lambda i: x[i]); r = [0.0] * len(x); i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and x[o[j + 1]] == x[o[i]]:
                j += 1
            for k in range(i, j + 1):
                r[o[k]] = (i + j) / 2
            i = j + 1
        return r
    rx, ry = rank(xs), rank(ys)
    n = len(rx); mx, my = sum(rx) / n, sum(ry) / n
    rho = sum((rx[i] - mx) * (ry[i] - my) for i in range(n)) / (
        sum((v - mx) ** 2 for v in rx) * sum((v - my) ** 2 for v in ry)) ** 0.5
    print(f"Spearman ΔΔG x dentes ausentes (n={n}): {rho:.2f}")


if __name__ == "__main__":
    main()
