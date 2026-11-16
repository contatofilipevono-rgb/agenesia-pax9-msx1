"""Consolida variantes missense de PAX9 e MSX1: ClinVar + gnomAD + AlphaMissense.

Uso: python -I consolidar.py <pasta_dados_brutos> <saida.csv>
Rótulos de referência (benchmark):
  P = Pathogenic / Likely pathogenic no ClinVar (sem conflito)
  B = Benign / Likely benign no ClinVar, ou frequência alélica no gnomAD >= 1e-3
  U = resto (incerto, conflitante ou só raro no gnomAD)
"""
import csv
import re
import sys
from pathlib import Path

DOMINIO = {"PAX9": (4, 130, "paired"), "MSX1": (172, 231, "homeobox")}
AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
       "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
       "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
MIS = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def ler(p):
    with open(p, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main():
    pasta, saida = Path(sys.argv[1]), Path(sys.argv[2])
    linhas = []
    for gene, (ini, fim, dom) in DOMINIO.items():
        am = {r["protein_variant"]: r for r in ler(pasta / f"alphamissense_{gene}.csv")}
        var = {}
        for r in ler(pasta / f"clinvar_{gene}.csv"):
            if "missense" not in r["molecular_consequence"]:
                continue
            for pc in r["protein_change"].split(","):
                pc = pc.strip()
                if MIS.match(pc) and pc in am:  # só a isoforma canônica casa com o AlphaMissense
                    v = var.setdefault(pc, {})
                    v.setdefault("clinvar", []).append(r)
                    break
        for r in ler(pasta / f"gnomad_{gene}.csv"):
            if r["consequence"] != "missense_variant":
                continue
            m = re.match(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", r["hgvsp"])
            if not m or m.group(1) not in AA3 or m.group(3) not in AA3:
                continue
            pc = AA3[m.group(1)] + m.group(2) + AA3[m.group(3)]
            if pc not in am:
                continue
            v = var.setdefault(pc, {})
            v["gnomad_af"] = max(float(r["af"] or 0), v.get("gnomad_af", 0))
            v["gnomad_ac"] = v.get("gnomad_ac", 0) + int(r["ac"])
        lit_path = pasta.parent / "literatura_missense.csv"
        lit = {r["variante"]: r for r in ler(lit_path) if r["gene"] == gene} if lit_path.exists() else {}
        for pc in lit:
            if pc in am:
                var.setdefault(pc, {})
        for pc, v in var.items():
            pos = int(MIS.match(pc).group(2))
            cls = [c["classification"] for c in v.get("clinvar", [])]
            rev = [c["review_status"] for c in v.get("clinvar", [])]
            cond = sorted({x for c in v.get("clinvar", []) for x in c["conditions"].split(";") if x})
            joined = " | ".join(cls)
            af = v.get("gnomad_af", 0)
            if re.search(r"(?<!Likely )Pathogenic|Likely pathogenic", joined) and "onflicting" not in joined \
                    and "enign" not in joined:
                rot = "P"
            elif ("enign" in joined and "athogenic" not in joined) or af >= 1e-3:
                rot = "B"
            elif pc in lit and "onflicting" not in joined and af < 1e-4:
                rot = "P-lit"  # descrita em família com agenesia, sem classificação P no ClinVar
            else:
                rot = "U"
            linhas.append({
                "gene": gene, "variante": pc, "posicao": pos,
                "no_dominio": f"sim ({dom})" if ini <= pos <= fim else "nao",
                "clinvar": joined, "clinvar_revisao": " | ".join(rev),
                "condicoes": "; ".join(cond),
                "gnomad_af": f"{af:.2e}" if af else "", "gnomad_ac": v.get("gnomad_ac", ""),
                "alphamissense": am[pc]["am_pathogenicity"], "am_classe": am[pc]["am_class"],
                "rotulo": rot,
                "literatura": (lit[pc]["referencia_original"] + " (via " + lit[pc]["fonte_revisao"] + ")")
                if pc in lit else "",
                "clinvar_ids": " ".join(c["accession"] for c in v.get("clinvar", [])),
            })
    linhas.sort(key=lambda r: (r["gene"], r["posicao"], r["variante"]))
    with open(saida, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(linhas[0].keys()))
        w.writeheader()
        w.writerows(linhas)

    # resumo
    for gene in DOMINIO:
        g = [r for r in linhas if r["gene"] == gene]
        for rot in ("P", "P-lit", "B", "U"):
            s = [r for r in g if r["rotulo"] == rot]
            d = [r for r in s if r["no_dominio"] != "nao"]
            print(f"{gene} {rot}: {len(s)} total, {len(d)} no domínio de ligação ao DNA")
    # desempenho do AlphaMissense (sem perda do códon de início)
    def avaliar(nome, positivos):
        ok = [r for r in linhas if not r["variante"].startswith("M1")]
        pos = [float(r["alphamissense"]) for r in ok if r["rotulo"] in positivos]
        neg = [float(r["alphamissense"]) for r in ok if r["rotulo"] == "B"]
        if pos and neg:
            auc = sum((p > n) + 0.5 * (p == n) for p in pos for n in neg) / (len(pos) * len(neg))
            tp = sum(s > 0.564 for s in pos)
            print(f"{nome}: AUC {auc:.3f} ({len(pos)} pos, {len(neg)} B); pos com am>0.564: {tp}/{len(pos)}")

    avaliar("ClinVar P", ("P",))
    avaliar("ClinVar P + literatura", ("P", "P-lit"))
    lit_only = [r for r in linhas if r["rotulo"] == "P-lit"]
    for r in lit_only:
        print("  P-lit", r["gene"], r["variante"], r["no_dominio"], r["alphamissense"], r["am_classe"])


if __name__ == "__main__":
    main()
