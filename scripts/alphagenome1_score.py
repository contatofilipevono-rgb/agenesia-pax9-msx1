"""Pontua as mesmas variantes no AlphaGenome 1 (API v1) para comparar com o AlphaGenome 2.

Uso: <venv-python> alphagenome1_score.py <alphagenome_todas.csv> <.env> <saida.csv>
O v1 não pede condição: devolve todos os tecidos/células. Resumimos de dois jeitos por tipo de saída:
  *_todos: maior |escore| em qualquer trilha do próprio gene;
  *_3cond: só nas mesmas 3 condições usadas no v2 (crista neural, mesenquimal, língua), se existirem no v1.
"""
import csv
import sys

import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v1 import dna_client
from alphagenome.models.v1 import variant_scorers as vs

TIPOS = ["SPLICE_SITES", "SPLICE_SITE_USAGE", "SPLICE_JUNCTIONS", "RNA_SEQ"]
CONDS = ("CL:0011012", "CL:0000134", "UBERON:0001723")


def chave(env):
    for ln in open(env, encoding="utf-8"):
        if ln.startswith("ALPHAGENOME_API_KEY="):
            return ln.split("=", 1)[1].strip()
    raise SystemExit("chave não encontrada")


def main():
    entrada, env, saida = sys.argv[1:4]
    with open(entrada, encoding="utf-8-sig") as f:
        rows = [r for r in csv.DictReader(f)]
    model = dna_client.create(chave(env))
    scorers = [vs.RECOMMENDED_VARIANT_SCORERS[t] for t in TIPOS]
    variantes = [genome.Variant(chromosome=r["cromossomo"], position=int(r["posicao_1based"]),
                                reference_bases=r["ref"], alternate_bases=r["alt"]) for r in rows]
    intervalos = [v.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB) for v in variantes]
    res = model.score_variants(intervalos, variantes, variant_scorers=scorers, max_workers=2)
    brutos, out = [], []
    for r, scores in zip(rows, res):
        df = vs.tidy_scores(scores)
        df["accession"] = r["accession"]
        brutos.append(df)
        if "gene_name" in df:
            df = df[df["gene_name"].isna() | (df["gene_name"] == r["gene"])]
        col = "quantile_score" if "quantile_score" in df else "raw_score"
        curie = next((c for c in ("ontology_curie", "biosample_curie") if c in df), None)
        linha = {k: r[k] for k in ("gene", "accession", "title", "grupo")}
        for t in TIPOS:
            sub = df[df["output_type"].astype(str).str.upper().str.endswith(t)]
            linha[f"v1_{t.lower()}_todos"] = round(float(sub[col].abs().max()), 4) if len(sub) else ""
            if curie:
                s3 = sub[sub[curie].astype(str).isin(CONDS)]
                linha[f"v1_{t.lower()}_3cond"] = round(float(s3[col].abs().max()), 4) if len(s3) else ""
        out.append(linha)
    pd.concat(brutos).to_csv(saida.replace(".csv", "_bruto.csv"), index=False, encoding="utf-8-sig")
    pd.DataFrame(out).to_csv(saida, index=False, encoding="utf-8-sig")
    print("coluna de escore:", col, "| coluna de célula:", curie, "| variantes:", len(out))


if __name__ == "__main__":
    main()
