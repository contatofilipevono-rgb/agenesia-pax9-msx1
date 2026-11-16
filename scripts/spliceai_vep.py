"""Busca escores SpliceAI (pré-calculados) pelo Ensembl VEP REST e compara com o AlphaGenome.

Uso: python -I spliceai_vep.py <alphagenome_todas.csv> <saida.csv>
SpliceAI = maior DS (AG, AL, DG, DL) entre os transcritos do próprio gene. Também registra as consequências do VEP.
"""
import csv
import json
import sys
import time
import urllib.request

URL = "https://rest.ensembl.org/vep/human/region?SpliceAI=1"


def vep(lines):
    body = json.dumps({"variants": lines}).encode()
    req = urllib.request.Request(URL, data=body, headers={
        "Content-Type": "application/json", "Accept": "application/json",
        "User-Agent": "academic-variant-curation/0.1"})
    for i in range(4):
        try:
            return json.loads(urllib.request.urlopen(req, timeout=180).read())
        except Exception:  # noqa: BLE001
            if i == 3:
                raise
            time.sleep(5 * (i + 1))


def auc(pos, neg):
    return sum((p > n) + 0.5 * (p == n) for p in pos for n in neg) / (len(pos) * len(neg))


def main():
    entrada, saida = sys.argv[1:3]
    with open(entrada, encoding="utf-8-sig") as f:
        rows = [r for r in csv.DictReader(f) if r["max_splice_sites_q"] not in ("", "nan")]
    lines, ids = [], {}
    for r in rows:
        chrom = r["cromossomo"].replace("chr", "")
        vid = f"v{len(lines)}"
        lines.append(f"{chrom} {r['posicao_1based']} {vid} {r['ref']} {r['alt']} . . .")
        ids[vid] = r
    out = vep(lines)
    for res in out:
        r = ids.get(res.get("id"))
        if not r:
            continue
        ds, cons = [], set()
        for tc in res.get("transcript_consequences", []):
            if tc.get("gene_symbol") != r["gene"]:
                continue
            cons.update(tc.get("consequence_terms", []))
            s = tc.get("spliceai")
            if s:
                ds.append(max(s.get(k, 0) for k in ("DS_AG", "DS_AL", "DS_DG", "DS_DL")))
        r["spliceai_max"] = max(ds) if ds else ""
        r["vep_consequencias"] = ";".join(sorted(cons))
    campos = list(rows[0].keys()) + [c for c in ("spliceai_max", "vep_consequencias") if c not in rows[0]]
    with open(saida, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=campos, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    # comparação (exclui perda de códon de início, que não é mecanismo de splicing)
    ok = [r for r in rows if r.get("spliceai_max") not in ("", None) and "start_lost" not in r.get("vep_consequencias", "")]
    for nome, key in (("AlphaGenome splicing", "max_splice_sites_q"), ("SpliceAI", "spliceai_max")):
        P = [float(r[key]) for r in ok if r["grupo"] == "P"]
        B = [float(r[key]) for r in ok if r["grupo"] == "B"]
        print(f"{nome}: AUC {auc(P, B):.3f} ({len(P)} P, {len(B)} B); max B {max(B):.3f}; min P {min(P):.3f}")
    sem = [r["title"][:40] for r in rows if r.get("spliceai_max") in ("", None)]
    print("sem SpliceAI:", len(sem), sem[:6])
    print("\nIncertas (U) ordenadas pelo AlphaGenome:")
    for r in sorted([r for r in ok if r["grupo"] == "U"], key=lambda r: -float(r["max_splice_sites_q"]))[:14]:
        print(f"  {r['gene']} {r['title'][:46]:46s} AG={float(r['max_splice_sites_q']):.3f} SpliceAI={float(r['spliceai_max']):.2f}")


if __name__ == "__main__":
    main()
