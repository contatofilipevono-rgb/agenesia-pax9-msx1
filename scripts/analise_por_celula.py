"""Análise por célula: mesenquimal como principal; crista neural, língua e osteoblasto como sensibilidade.

Uso: python -I analise_por_celula.py <splicing_alphagenome_vs_spliceai.csv> <alphagenome_osteoblasto.csv> <saida.csv>
Exclui perda do códon de início. Limiar de cada célula = maior escore entre as benignas naquela célula.
"""
import csv
import sys

CEL = ["mesenquimal", "crista_neural", "lingua", "osteoblasto"]
TIPOS = ["splice_sites", "rna_seq"]


def f(x):
    try:
        v = float(x)
        return v if v == v else None
    except (TypeError, ValueError):
        return None


def auc(P, B):
    return sum((p > b) + 0.5 * (p == b) for p in P for b in B) / (len(P) * len(B))


def main():
    a, ost, saida = sys.argv[1:4]
    rows = [r for r in csv.DictReader(open(a, encoding="utf-8-sig"))
            if "start_lost" not in r.get("vep_consequencias", "")]
    o = {r["accession"]: r for r in csv.DictReader(open(ost, encoding="utf-8-sig"))}
    for r in rows:
        for t in TIPOS:
            r[f"osteoblasto_{t}_q"] = o.get(r["accession"], {}).get(f"osteoblasto_{t}_q", "")
    print(f"{'célula':14s} {'AUC splicing':>13s} {'AUC expressão':>14s} {'limiar splicing (B máx)':>24s} incertas acima")
    for c in CEL:
        lin = []
        for t in TIPOS:
            k = f"{c}_{t}_q"
            P = [f(r[k]) for r in rows if r["grupo"] == "P" and f(r[k]) is not None]
            B = [f(r[k]) for r in rows if r["grupo"] == "B" and f(r[k]) is not None]
            lin.append((auc(P, B), max(B)))
        k = f"{c}_splice_sites_q"
        acima = [r for r in rows if r["grupo"] == "U" and f(r[k]) is not None and f(r[k]) > lin[0][1]]
        for r in rows:
            r[f"{c}_acima_limiar"] = "sim" if r in acima else ""
        print(f"{c:14s} {lin[0][0]:13.3f} {lin[1][0]:14.3f} {lin[0][1]:24.3f} {len(acima)}: "
              + ", ".join(r["title"].split(":")[-1][:14] for r in acima))
    with open(saida, "w", newline="", encoding="utf-8-sig") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    main()
