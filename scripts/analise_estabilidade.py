"""Cruza o ΔΔG do ThermoMPNN com as mutações patogênicas e a posição em relação ao DNA.

Uso: python -I analise_estabilidade.py <pasta_estabilidade> <resultados_afserver.csv>
- Convenção: o ThermoMPNN não deixa explícito no script; inferimos pelo sinal da mediana de todas as trocas
  (a maioria das trocas ao acaso desestabiliza). Reportamos 'desestab' = ΔΔG com sinal tal que >0 desestabiliza.
- Percentil de cada mutação entre todas as trocas possíveis do domínio.
- Compara: patogênicas x todas as trocas possíveis (Mann-Whitney por U/AUC); contato x distantes.
"""
import csv
import os
import sys
from statistics import median

INI = {"PAX9": 6, "MSX1": 166}


def auc(a, b):
    return sum((x > y) + 0.5 * (x == y) for x in a for y in b) / (len(a) * len(b))


def main():
    pasta, res_csv = sys.argv[1:3]
    ssm, wt_aa = {}, {}
    for gene in INI:
        path = os.path.join(pasta, f"ThermoMPNN_inference_{gene}_WT_cadeiaA.csv")
        for r in csv.DictReader(open(path, encoding="utf-8")):
            if r["wildtype"] == r["mutation"]:
                continue
            # ATENÇÃO: o ThermoMPNN numera posições a partir de 0
            ssm[(gene, int(r["position"]), r["mutation"])] = float(r["ddG_pred"])
            wt_aa[(gene, int(r["position"]))] = r["wildtype"]
    sinal = {}
    for gene in INI:
        vals = [v for (g, _, _), v in ssm.items() if g == gene]
        med = median(vals)
        sinal[gene] = 1 if med > 0 else -1  # multiplica para que >0 = desestabilizante
        print(f"{gene}: {len(vals)} trocas; mediana bruta {med:.2f} -> sinal {'+' if sinal[gene] > 0 else '-'} desestabiliza")
    rows = list(csv.DictReader(open(res_csv, encoding="utf-8-sig")))
    for r in rows:
        r["ddg_thermompnn"] = r["desestab_pct"] = ""
        if r["variante"] == "WT":
            continue
        g, v = r["gene"], r["variante"]
        pos = int(v[1:-1]) - INI[g]  # 0-based na construção
        assert wt_aa.get((g, pos)) == v[0], (g, v, wt_aa.get((g, pos)))
        d = ssm.get((g, pos, v[-1]))
        if d is None:
            continue
        d *= sinal[g]
        todos = [x * sinal[g] for (gg, _, _), x in ssm.items() if gg == g]
        r["ddg_thermompnn"] = round(d, 2)
        r["desestab_pct"] = round(100 * sum(t < d for t in todos) / len(todos), 1)
    with open(res_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)

    for g in INI:
        mut = [r for r in rows if r["gene"] == g and r["ddg_thermompnn"] != ""]
        todos = [x * sinal[g] for (gg, _, _), x in ssm.items() if gg == g]
        pat = [float(r["ddg_thermompnn"]) for r in mut]
        print(f"\n{g}: {len(pat)} patogênicas | mediana ΔΔG {median(pat):.2f} vs todas as trocas {median(todos):.2f}"
              f" | AUC patogênicas > todas = {auc(pat, todos):.2f}")
        for classe in ("contato DNA", "próximo", "distante"):
            c = [float(r["ddg_thermompnn"]) for r in mut if r["classe_sitio"] == classe]
            if c:
                print(f"   {classe:12s} n={len(c):2d} mediana ΔΔG {median(c):.2f}")
        for r in sorted(mut, key=lambda r: -float(r["ddg_thermompnn"])):
            print(f"   {r['variante']:6s} {r['classe_sitio']:12s} ΔΔG={float(r['ddg_thermompnn']):5.2f} "
                  f"(percentil {r['desestab_pct']}) dist={r['dist_dna_wt']}")


if __name__ == "__main__":
    main()
