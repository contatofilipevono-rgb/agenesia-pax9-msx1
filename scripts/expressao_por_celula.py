"""Checa se PAX9 e MSX1 são expressos (previsão de RNA-seq do AlphaGenome 2) em cada célula usada.

Para cada célula e gene: prevê RNA-seq numa janela de 100 kb centrada no gene e tira a média do sinal
dentro do gene (maior média entre as trilhas). Controles: GAPDH (expresso em quase toda célula) e ALB (fígado).
Uso: <venv-python> expressao_por_celula.py <.env> <saida.csv>
"""
import sys

import numpy as np
import pandas as pd
from alphagenome.data import genome
from alphagenome.models.v2 import dna_client

GENES = {  # GRCh38, Ensembl REST (15/11/2026)
    "PAX9": ("chr14", 36657568, 36679362),
    "MSX1": ("chr4", 4859665, 4863936),
    "GAPDH": ("chr12", 6534012, 6538930),
    "ALB": ("chr4", 73397114, 73426607),
}
CELULAS = {"CL:0000134": "mesenquimal", "CL:0011012": "crista_neural", "UBERON:0001723": "lingua",
           "UBERON:0002481_CL:0000062": "osteoblasto"}


def chave(env):
    for ln in open(env, encoding="utf-8"):
        if ln.startswith("ALPHAGENOME_API_KEY="):
            return ln.split("=", 1)[1].strip()


def main():
    env, saida = sys.argv[1:3]
    model = dna_client.create(chave(env))
    linhas = []
    for cond, nome in CELULAS.items():
        for g, (ch, ini, fim) in GENES.items():
            meio = (ini + fim) // 2
            iv = genome.Interval(ch, meio - dna_client.SEQUENCE_LENGTH_100KB // 2,
                                 meio + dna_client.SEQUENCE_LENGTH_100KB // 2)
            out = model.predict_interval(iv, requested_outputs=[dna_client.OutputType.RNA_SEQ], condition=cond)
            td = out.rna_seq
            vals = np.asarray(td.values)  # posições x trilhas
            res = vals.shape[0] / iv.width  # posições por pb
            a = int((ini - iv.start) * res)
            b = int((fim - iv.start) * res)
            media = float(vals[a:b].mean(axis=0).max()) if b > a else float("nan")
            pico = float(np.percentile(vals[a:b], 99, axis=0).max()) if b > a else float("nan")  # sinal dos éxons
            linhas.append({"celula": nome, "condicao": cond, "gene": g, "rna_media_no_gene": round(media, 4),
                           "rna_p99_no_gene": round(pico, 4), "n_trilhas": vals.shape[1]})
            print(nome, g, round(media, 4), round(pico, 4))
    df = pd.DataFrame(linhas)
    # normaliza pelo GAPDH da mesma célula
    gap = df[df.gene == "GAPDH"].set_index("celula")["rna_p99_no_gene"]
    alb = df[df.gene == "ALB"].set_index("celula")["rna_p99_no_gene"]
    df["p99_relativo_gapdh"] = (df["rna_p99_no_gene"] / df["celula"].map(gap)).round(4)
    df["p99_vezes_alb"] = (df["rna_p99_no_gene"] / df["celula"].map(alb)).round(1)
    df.to_csv(saida, index=False, encoding="utf-8-sig")
    print(df.pivot(index="gene", columns="celula", values="rna_p99_no_gene").to_string())
    print(df.pivot(index="gene", columns="celula", values="p99_vezes_alb").to_string())


if __name__ == "__main__":
    main()
