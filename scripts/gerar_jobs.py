"""Gera as sequências dos jobs do AlphaFold Server (proteína normal e mutantes + DNA).

Construções:
  PAX9 resíduos 6–137 (domínio paired) + DNA da estrutura 6PAX (PAX6, parente).
  MSX1 resíduos 166–237 (homeobox 172–231 com flancos) + DNA da estrutura 1IG7 (Msx-1 de camundongo).
Uso: python -I gerar_jobs.py <variantes.csv> <saida.csv>
"""
import csv
import sys
import urllib.request

CONSTR = {
    "PAX9": {"uniprot": "P55771", "ini": 6, "fim": 137,
             "dna": ["AAGCATTTTCACGCATGAGTGCACAG", "TTCTGTGCACTCATGCGTGAAAATGC"]},
    "MSX1": {"uniprot": "P28360", "ini": 166, "fim": 237,
             "dna": ["TCCTTCAATTAGTGAC", "CACTAATTGAAGG"]},
}


def seq_uniprot(acc):
    req = urllib.request.Request(f"https://rest.uniprot.org/uniprotkb/{acc}.fasta",
                                 headers={"User-Agent": "academic-variant-curation/0.1"})
    txt = urllib.request.urlopen(req, timeout=60).read().decode()
    return "".join(l.strip() for l in txt.splitlines() if not l.startswith(">"))


def main():
    variantes, saida = sys.argv[1], sys.argv[2]
    with open(variantes, encoding="utf-8-sig") as f:
        linhas = list(csv.DictReader(f))
    jobs = []
    for gene, c in CONSTR.items():
        full = seq_uniprot(c["uniprot"])
        wt = full[c["ini"] - 1:c["fim"]]
        jobs.append({"job": f"{gene}_WT", "gene": gene, "variante": "WT", "rotulo": "",
                     "proteina": wt, "dna1": c["dna"][0], "dna2": c["dna"][1]})
        sel = [r for r in linhas if r["gene"] == gene and r["rotulo"] in ("P", "P-lit")
               and r["no_dominio"] != "nao"]
        sel.sort(key=lambda r: int(r["posicao"]))
        for r in sel:
            v = r["variante"]
            ref, pos, alt = v[0], int(v[1:-1]), v[-1]
            assert full[pos - 1] == ref, (gene, v, full[pos - 1])
            assert c["ini"] <= pos <= c["fim"], (gene, v)
            i = pos - c["ini"]
            mut = wt[:i] + alt + wt[i + 1:]
            jobs.append({"job": f"{gene}_{v}", "gene": gene, "variante": v, "rotulo": r["rotulo"],
                         "proteina": mut, "dna1": c["dna"][0], "dna2": c["dna"][1]})
    with open(saida, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(jobs[0].keys()) + ["status", "data_envio", "iptm", "ptm"])
        w.writeheader()
        for j in jobs:
            w.writerow(j)
    for j in jobs:
        print(j["job"], len(j["proteina"]))
    print("total", len(jobs))


if __name__ == "__main__":
    main()
