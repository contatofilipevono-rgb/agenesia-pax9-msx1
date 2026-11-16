"""Lê o .zip em lote do AlphaFold Server (sem extrair nem executar nada) e resume cada job.

Para cada job escolhe, entre as 5 amostras, a de maior ranking_score e registra:
  iptm, ptm, ranking_score, has_clash, fraction_disordered,
  ipTM proteína–DNA (média dos pares A-B e A-C), PAE mínimo proteína–DNA,
  pLDDT médio do resíduo mutado e probabilidade máxima de contato desse resíduo com o DNA.
Uso: python -I extrair_afserver.py <lote.zip> <jobs_afserver.csv> <saida.csv>
"""
import csv
import json
import re
import sys
import zipfile
from statistics import mean

CONSTR_INI = {"PAX9": 6, "MSX1": 166}


def main():
    zpath, jobs_csv, saida = sys.argv[1:4]
    with open(jobs_csv, encoding="utf-8-sig") as f:
        jobs = {r["job"].lower(): r for r in csv.DictReader(f)}
    jobs["pax9_paired_domain_dna_wt"] = jobs.get("pax9_wt")
    z = zipfile.ZipFile(zpath)
    nomes = z.namelist()
    pastas = sorted({n.split("/")[0] for n in nomes if "/" in n})
    linhas = []
    for p in pastas:
        job = jobs.get(p)
        if not job:
            print("sem job para", p)
            continue
        melhores = []
        for k in range(5):
            s = f"{p}/fold_{p}_summary_confidences_{k}.json"
            if s in nomes:
                melhores.append((json.loads(z.read(s))["ranking_score"], k))
        rs, k = max(melhores)
        summ = json.loads(z.read(f"{p}/fold_{p}_summary_confidences_{k}.json"))
        full = json.loads(z.read(f"{p}/fold_{p}_full_data_{k}.json"))
        cpi = summ["chain_pair_iptm"]
        cpp = summ["chain_pair_pae_min"]
        gene, var = job["gene"], job["variante"]
        n_prot = len(job["proteina"])
        tok_chain = full["token_chain_ids"]
        prot_tokens = [i for i, c in enumerate(tok_chain) if c == "A"]
        dna_tokens = [i for i, c in enumerate(tok_chain) if c != "A"]
        assert len(prot_tokens) == n_prot, (p, len(prot_tokens), n_prot)
        site_plddt = site_contact = ""
        if var != "WT":
            pos = int(re.match(r"[A-Z](\d+)[A-Z]", var).group(1))
            idx = prot_tokens[pos - CONSTR_INI[gene]]
            cp = full["contact_probs"][idx]
            site_contact = round(max(cp[j] for j in dna_tokens), 3)
            # pLDDT do resíduo: média do B-factor (pLDDT) dos átomos desse resíduo no mmCIF do modelo
            cif = z.read(f"{p}/fold_{p}_model_{k}.cif").decode()
            cols, pls = [], []
            resnum = pos - CONSTR_INI[gene] + 1  # numeração 1-based dentro da construção
            for ln in cif.splitlines():
                if ln.startswith("_atom_site."):
                    cols.append(ln.split(".", 1)[1].strip())
                elif ln.startswith(("ATOM", "HETATM")) and cols:
                    f = ln.split()
                    rec = dict(zip(cols, f))
                    if rec.get("label_asym_id") == "A" and rec.get("label_seq_id") == str(resnum):
                        pls.append(float(rec["B_iso_or_equiv"]))
            if pls:
                site_plddt = round(mean(pls), 1)
        linhas.append({
            "job": job["job"], "gene": gene, "variante": var, "rotulo": job["rotulo"], "amostra": k,
            "ranking_score": round(rs, 3), "iptm": summ["iptm"], "ptm": summ["ptm"],
            "iptm_prot_dna": round(mean([cpi[0][1], cpi[0][2]]), 3),
            "pae_min_prot_dna": round(min(cpp[0][1], cpp[0][2]), 2),
            "has_clash": summ["has_clash"], "fraction_disordered": summ["fraction_disordered"],
            "plddt_sitio": site_plddt, "contato_sitio_dna": site_contact,
        })
    linhas.sort(key=lambda r: (r["gene"], r["variante"] != "WT", r["job"]))
    with open(saida, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()))
        w.writeheader()
        w.writerows(linhas)
    for r in linhas:
        print(r["job"], r["iptm"], r["ptm"], r["iptm_prot_dna"], r["pae_min_prot_dna"], r["contato_sitio_dna"],
              r["plddt_sitio"], r["has_clash"])


if __name__ == "__main__":
    main()
