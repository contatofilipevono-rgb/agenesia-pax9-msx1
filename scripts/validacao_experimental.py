"""Compara as previsões (contato com DNA, ΔΔG ThermoMPNN) com os experimentos de Wang et al. 2009 (Hum Mol Genet,
doi:10.1093/hmg/ddp221; PMC2706687): EMSA, estabilidade modelada pelos autores e gravidade (dentes ausentes).
Uso: python -I -X utf8 validacao_experimental.py <validacao_wang2009.csv> <resultados_afserver.csv> [corte_ddg=1.5]
"""
import csv
import sys


def rank(x):
    o = sorted(range(len(x)), key=lambda i: x[i])
    r = [0.0] * len(x)
    i = 0
    while i < len(o):
        j = i
        while j + 1 < len(o) and x[o[j + 1]] == x[o[i]]:
            j += 1
        for k in range(i, j + 1):
            r[o[k]] = (i + j) / 2
        i = j + 1
    return r


def spearman(a, b):
    ra, rb = rank(a), rank(b)
    n = len(a)
    ma, mb = sum(ra) / n, sum(rb) / n
    num = sum((ra[i] - ma) * (rb[i] - mb) for i in range(n))
    return num / (sum((v - ma) ** 2 for v in ra) * sum((v - mb) ** 2 for v in rb)) ** 0.5


def main():
    exp = {r["variante"]: r for r in csv.DictReader(open(sys.argv[1], encoding="utf-8"))}
    res = {r["variante"]: r for r in csv.DictReader(open(sys.argv[2], encoding="utf-8-sig")) if r["gene"] == "PAX9"}
    vs = list(exp)
    corte = float(sys.argv[3]) if len(sys.argv) > 3 else 1.5  # 16/11: padrão passou de 1 para 1,5 (corte do artigo)
    print(f"regra: contato < 4 Å ou ΔΔG >= {corte}")
    print(f"{'var':5s} {'posição':12s} {'ΔΔG':>5s} {'previsto':>9s} {'EMSA':>5s} {'observado':>9s} dentes")
    ac = 0
    for v in vs:
        r, e = res[v], exp[v]
        ddg = float(r["ddg_thermompnn"])
        prev = "perde" if (r["classe_sitio"] == "contato DNA" or ddg >= corte) else "mantém"
        obs = "perde" if float(e["ligacao_escore"]) <= 1 else "mantém"
        ac += prev == obs
        print(f"{v:5s} {r['classe_sitio']:12s} {ddg:5.2f} {prev:>9s} {e['ligacao_dna_emsa']:>5s} {obs:>9s} "
              f"{e['dentes_ausentes_medio']}")
    print(f"acertos na ligação ao DNA: {ac}/{len(vs)}")
    ddg = [float(res[v]["ddg_thermompnn"]) for v in vs]
    wang = [-float(exp[v]["estabilidade_modelo_wang"]) for v in vs]
    lig = [float(exp[v]["ligacao_escore"]) for v in vs]
    dentes = [float(exp[v]["dentes_ausentes_medio"]) for v in vs]
    print(f"Spearman ΔΔG ThermoMPNN x desestabilização modelada por Wang: {spearman(ddg, wang):.2f}")
    print(f"Spearman ΔΔG x ligação ao DNA (EMSA): {spearman(ddg, lig):.2f}")
    print(f"Spearman ligação ao DNA x dentes ausentes: {spearman(lig, dentes):.2f}")
    print(f"Spearman ΔΔG x dentes ausentes: {spearman(ddg, dentes):.2f}")


if __name__ == "__main__":
    main()
