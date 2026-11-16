"""Compara AlphaGenome 1, AlphaGenome 2 e SpliceAI nas mesmas variantes não codificantes.

Uso: python -I comparar_versoes.py <splicing_alphagenome_vs_spliceai.csv> <alphagenome1_todas.csv> <saida.csv>
"""
import csv
import sys


def f(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def auc(rows, key):
    P = [f(r[key]) for r in rows if r["grupo"] == "P" and f(r[key]) is not None]
    B = [f(r[key]) for r in rows if r["grupo"] == "B" and f(r[key]) is not None]
    if not P or not B:
        return None, len(P), len(B)
    return sum((p > b) + 0.5 * (p == b) for p in P for b in B) / (len(P) * len(B)), len(P), len(B)


def spearman(rows, a, b):
    pares = [(f(r[a]), f(r[b])) for r in rows if f(r[a]) is not None and f(r[b]) is not None]

    def rank(xs):
        o = sorted(range(len(xs)), key=lambda i: xs[i])
        rk = [0.0] * len(xs)
        i = 0
        while i < len(o):
            j = i
            while j + 1 < len(o) and xs[o[j + 1]] == xs[o[i]]:
                j += 1
            for k in range(i, j + 1):
                rk[o[k]] = (i + j) / 2
            i = j + 1
        return rk
    x, y = rank([p[0] for p in pares]), rank([p[1] for p in pares])
    n = len(x)
    mx, my = sum(x) / n, sum(y) / n
    cov = sum((x[i] - mx) * (y[i] - my) for i in range(n))
    return cov / ((sum((v - mx) ** 2 for v in x) * sum((v - my) ** 2 for v in y)) ** 0.5), n


def main():
    a, b, saida = sys.argv[1:4]
    v2 = {r["accession"]: r for r in csv.DictReader(open(a, encoding="utf-8-sig"))}
    rows = []
    for r in csv.DictReader(open(b, encoding="utf-8-sig")):
        base = v2.get(r["accession"])
        if not base or "start_lost" in base.get("vep_consequencias", ""):
            continue
        base.update({k: v for k, v in r.items() if k.startswith("v1_")})
        rows.append(base)
    chaves = [("AG2 splicing (3 células)", "max_splice_sites_q"),
              ("AG1 splicing (todas)", "v1_splice_sites_todos"),
              ("AG1 splicing (3 células)", "v1_splice_sites_3cond"),
              ("SpliceAI", "spliceai_max"),
              ("AG2 expressão (3 células)", "max_rna_seq_q"),
              ("AG1 expressão (todas)", "v1_rna_seq_todos"),
              ("AG1 expressão (3 células)", "v1_rna_seq_3cond")]
    print("AUC patogênicas x benignas:")
    for nome, k in chaves:
        v, p, n = auc(rows, k)
        print(f"  {nome:28s} {'—' if v is None else f'{v:.3f}'}  ({p} P, {n} B)")
    print("Correlação de Spearman (todas as variantes):")
    for x, y in (("max_splice_sites_q", "v1_splice_sites_todos"), ("max_splice_sites_q", "spliceai_max"),
                 ("v1_splice_sites_todos", "spliceai_max"), ("max_rna_seq_q", "v1_rna_seq_todos")):
        rho, n = spearman(rows, x, y)
        print(f"  {x} x {y}: rho={rho:.2f} (n={n})")
    print("Incertas com maior splicing em qualquer método:")
    U = [r for r in rows if r["grupo"] == "U"]
    U.sort(key=lambda r: -(f(r["max_splice_sites_q"]) or 0))
    for r in U[:10]:
        print(f"  {r['gene']} {r['title'][:44]:44s} AG2={f(r['max_splice_sites_q']):.3f} "
              f"AG1={f(r['v1_splice_sites_todos']) or 0:.3f} SpliceAI={f(r['spliceai_max']) or 0:.2f}")
    with open(saida, "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
