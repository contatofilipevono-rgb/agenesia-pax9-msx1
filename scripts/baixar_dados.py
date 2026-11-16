"""Baixa dados públicos de variantes de PAX9 e MSX1.

Fontes: ClinVar (NCBI E-utilities), gnomAD v4 (API GraphQL), AlphaMissense (AlphaFold DB).
Uso: python -I baixar_dados.py <pasta_saida>
"""
import csv
import json
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path

GENES = {"PAX9": "P55771", "MSX1": "P28360"}
UA = {"User-Agent": "academic-variant-curation/0.1"}


def get(url, data=None, headers=None, tries=4):
    h = dict(UA)
    if headers:
        h.update(headers)
    for i in range(tries):
        try:
            req = urllib.request.Request(url, data=data, headers=h)
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            if i == tries - 1:
                raise
            time.sleep(2 * (i + 1))


def clinvar(gene, out):
    base = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
    term = f"{gene}[gene] AND single_gene[prop]"
    q = urllib.parse.urlencode({"db": "clinvar", "term": term, "retmax": 5000, "retmode": "json"})
    ids = json.loads(get(base + "esearch.fcgi?" + q))["esearchresult"]["idlist"]
    rows = []
    for i in range(0, len(ids), 200):
        chunk = ids[i:i + 200]
        q = urllib.parse.urlencode({"db": "clinvar", "id": ",".join(chunk), "retmode": "json"})
        res = json.loads(get(base + "esummary.fcgi?" + q))["result"]
        for uid in res.get("uids", []):
            d = res[uid]
            germ = d.get("germline_classification") or d.get("clinical_significance") or {}
            traits = germ.get("trait_set") or d.get("trait_set") or []
            rows.append({
                "gene": gene,
                "clinvar_uid": uid,
                "accession": d.get("accession", ""),
                "title": d.get("title", ""),
                "protein_change": d.get("protein_change", ""),
                "molecular_consequence": ";".join(d.get("molecular_consequence_list", []) or []),
                "variant_type": d.get("obj_type", ""),
                "classification": germ.get("description", ""),
                "review_status": germ.get("review_status", ""),
                "conditions": ";".join(t.get("trait_name", "") for t in traits),
            })
        time.sleep(0.4)
    with open(out / f"clinvar_{gene}.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["gene"])
        w.writeheader()
        w.writerows(rows)
    return len(rows)


GNOMAD_Q = """
query($g: String!) {
  gene(gene_symbol: $g, reference_genome: GRCh38) {
    canonical_transcript_id
    variants(dataset: gnomad_r4) {
      variant_id consequence hgvsp hgvsc transcript_id
      exome { ac an }
      genome { ac an }
    }
  }
}"""


def gnomad(gene, out):
    body = json.dumps({"query": GNOMAD_Q, "variables": {"g": gene}}).encode()
    res = json.loads(get("https://gnomad.broadinstitute.org/api", data=body,
                         headers={"Content-Type": "application/json"}))
    g = res["data"]["gene"]
    rows = []
    for v in g["variants"]:
        ac = sum((v.get(k) or {}).get("ac", 0) or 0 for k in ("exome", "genome"))
        an = sum((v.get(k) or {}).get("an", 0) or 0 for k in ("exome", "genome"))
        rows.append({
            "gene": gene, "variant_id": v["variant_id"], "consequence": v["consequence"],
            "hgvsp": v.get("hgvsp") or "", "hgvsc": v.get("hgvsc") or "",
            "transcript": v.get("transcript_id") or "", "ac": ac, "an": an,
            "af": (ac / an) if an else "",
        })
    with open(out / f"gnomad_{gene}.csv", "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    return len(rows), g["canonical_transcript_id"]


def alphamissense(gene, uniprot, out):
    meta = json.loads(get(f"https://alphafold.ebi.ac.uk/api/prediction/{uniprot}"))
    url = meta[0].get("amAnnotationsUrl")
    if not url:
        return 0
    data = get(url)
    (out / f"alphamissense_{gene}.csv").write_bytes(data)
    return data.count(b"\n") - 1


def main():
    out = Path(sys.argv[1])
    out.mkdir(parents=True, exist_ok=True)
    for gene, up in GENES.items():
        print(gene, "ClinVar:", clinvar(gene, out))
        print(gene, "gnomAD:", gnomad(gene, out))
        print(gene, "AlphaMissense:", alphamissense(gene, up, out))


if __name__ == "__main__":
    main()
