"""Estabilidade (ΔΔG) dos mutantes de PAX9 e MSX1 com o DDMut (servidor gratuito, API pública).

1) Extrai a cadeia proteica (A) do modelo WT de melhor ranking no .zip do AlphaFold Server e grava em PDB.
2) Converte as variantes para a numeração do modelo (PAX9: pos-5; MSX1: pos-165) e confere o resíduo de referência.
3) Envia a lista ao endpoint prediction_list e espera o resultado.
Uso: python -I ddmut_estabilidade.py <lote.zip> <resultados_afserver.csv> <pasta_saida>
"""
import csv
import json
import os
import subprocess
import sys
import time
import zipfile

API = "https://biosig.lab.uq.edu.au/ddmut/api/prediction_list"
INI = {"PAX9": 6, "MSX1": 166}
WT_DIR = {"PAX9": "pax9_paired_domain_dna_wt", "MSX1": "msx1_wt"}
TRES = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
        "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
        "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V"}


def pdb_da_cadeia_a(z, pasta, destino):
    melhores = []
    for k in range(5):
        melhores.append((json.loads(z.read(f"{pasta}/fold_{pasta}_summary_confidences_{k}.json"))["ranking_score"], k))
    k = max(melhores)[1]
    cif = z.read(f"{pasta}/fold_{pasta}_model_{k}.cif").decode()
    cols, linhas, seq = [], [], {}
    n = 0
    for ln in cif.splitlines():
        if ln.startswith("_atom_site."):
            cols.append(ln.split(".", 1)[1].strip())
        elif ln.startswith("ATOM") and cols:
            r = dict(zip(cols, ln.split()))
            if r["label_asym_id"] != "A":
                continue
            n += 1
            nome = r["label_atom_id"].strip('"')
            res = r["label_comp_id"]
            num = int(r["label_seq_id"])
            seq[num] = TRES.get(res, "X")
            el = r["type_symbol"]
            nome_fmt = f" {nome:<3s}" if len(nome) < 4 else nome
            linhas.append(f"ATOM  {n:5d} {nome_fmt} {res:3s} A{num:4d}    "
                          f"{float(r['Cartn_x']):8.3f}{float(r['Cartn_y']):8.3f}{float(r['Cartn_z']):8.3f}"
                          f"  1.00{float(r['B_iso_or_equiv']):6.2f}          {el:>2s}")
    with open(destino, "w") as f:
        f.write("\n".join(linhas) + "\nTER\nEND\n")
    return seq


def curl(*args):
    return subprocess.run(["curl.exe", "-s", *args], capture_output=True, text=True, timeout=300).stdout


def main():
    zpath, res_csv, out = sys.argv[1:4]
    os.makedirs(out, exist_ok=True)
    z = zipfile.ZipFile(zpath)
    rows = [r for r in csv.DictReader(open(res_csv, encoding="utf-8-sig")) if r["variante"] != "WT"]
    resultado = {}
    for gene in ("PAX9", "MSX1"):
        pdb = os.path.join(out, f"{gene}_WT_cadeiaA.pdb")
        seq = pdb_da_cadeia_a(z, WT_DIR[gene], pdb)
        muts = {}
        for r in rows:
            if r["gene"] != gene:
                continue
            v = r["variante"]
            ref, pos, alt = v[0], int(v[1:-1]), v[-1]
            num = pos - INI[gene] + 1
            assert seq[num] == ref, (gene, v, seq[num])
            muts[f"{ref}{num}{alt}"] = v
        # o endpoint de lista recusa o envio; usamos o de mutação única, uma por vez
        single = API.replace("prediction_list", "prediction_single")
        jobs = {}
        for m, v in muts.items():
            resp = curl(single, "-X", "POST", "-F", f"pdb_file=@{pdb}", "-F", f"mutation={m}", "-F", "chain=A")
            jobs[m] = json.loads(resp)["job_id"]
            time.sleep(3)
        print(gene, "enviados:", len(jobs))
        brutos = {}
        pend = dict(jobs)
        for _ in range(160):
            if not pend:
                break
            time.sleep(15)
            for m, job in list(pend.items()):
                r = curl("-G", single, "--data-urlencode", f"job_id={job}")
                if not r.strip() or '"RUNNING"' in r or "QUEUED" in r.upper():
                    continue
                d = json.loads(r)
                brutos[m] = d
                if "prediction" in d:
                    resultado[(gene, muts[m])] = d["prediction"]
                del pend[m]
        if pend:
            print(gene, "sem resultado:", list(pend))
        open(os.path.join(out, f"{gene}_ddmut.json"), "w").write(json.dumps(brutos, indent=1))
    for r in rows:
        r["ddg_ddmut"] = resultado.get((r["gene"], r["variante"]), "")
    todas = list(csv.DictReader(open(res_csv, encoding="utf-8-sig")))
    for t in todas:
        t["ddg_ddmut"] = resultado.get((t["gene"], t["variante"]), "")
    with open(res_csv, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(todas[0].keys()))
        w.writeheader()
        w.writerows(todas)
    for r in sorted(rows, key=lambda r: float(r["ddg_ddmut"] or 0)):
        print(r["job"], r.get("classe_sitio", ""), r["ddg_ddmut"])


if __name__ == "__main__":
    main()
