"""Melhoria 4: aplica a regra (contato com DNA OU ΔΔG>=1) às variantes incertas (U) dentro dos domínios.

Calibração pelas patogênicas x controle gnomAD recorrente (AC>=2):
  sensibilidade = P(regra+ | patogênica); falso-positivo = P(regra+ | controle)
  LR+ = sens / FP ; LR- = (1-sens) / (1-FP)
Saída: lista de incertas com regra, ΔΔG, distância, AlphaMissense e LR aplicável.
Uso: python -I -X utf8 classificar_incertas.py <variantes_PAX9_MSX1.csv> <grupo_controle.csv> <lote.zip> <pasta_estabilidade> <saida.csv> [corte_ddg=1.5]
(16/11: limiar padrão passou de 1 para 1.5, o validado por leave-one-out; a versão com 1 está em incertas_classificadas_ddg1.csv)
"""
import csv
import math
import os
import runpy
import sys
import zipfile

INI = {"PAX9": 6, "MSX1": 166}
DOM = {"PAX9": (6, 130), "MSX1": (172, 231)}


def main():
    var_csv, ctrl_csv, zpath, pasta, saida = sys.argv[1:6]
    corte = float(sys.argv[6]) if len(sys.argv) > 6 else 1.5  # limiar validado por leave-one-out
    print(f"regra: contato < 4 Å ou ΔΔG >= {corte}")
    ctrl = list(csv.DictReader(open(ctrl_csv, encoding="utf-8-sig")))
    calib = {}
    for g in INI:
        P = [r for r in ctrl if r["gene"] == g and r["grupo"] == "patogenica"]
        C = [r for r in ctrl if r["gene"] == g and r["grupo"] == "controle_gnomad" and int(r["gnomad_ac"]) >= 2]
        pos = lambda r: r["contato_dna"] == "sim" or (r["ddg_thermompnn"] != "" and float(r["ddg_thermompnn"]) >= corte)
        sens = sum(map(pos, P)) / len(P)
        fp = sum(map(pos, C)) / len(C)
        calib[g] = (sens, fp, sens / fp, (1 - sens) / (1 - fp))
        print(f"{g}: sensibilidade {sens:.2f}, falso-positivo {fp:.2f}, LR+ {sens/fp:.1f}, LR- {(1-sens)/(1-fp):.2f}"
              f"  (n={len(P)} P, {len(C)} controle)")
    dm = runpy.run_path(os.path.join(os.path.dirname(__file__), "distancia_dna.py"))
    z = zipfile.ZipFile(zpath)
    dist = {}
    for g, p in dm["WT_DIR"].items():
        prot, dna = dm["ler_modelo"](z, p)
        for i, at in prot.items():
            dist[(g, i + INI[g] - 1)] = min(math.dist(a, b) for a, _ in at for b in dna)
    ddg = {}
    for g in INI:
        for r in csv.DictReader(open(os.path.join(pasta, f"ThermoMPNN_inference_{g}_WT_cadeiaA.csv"), encoding="utf-8")):
            if r["wildtype"] != r["mutation"]:
                ddg[(g, int(r["position"]) + INI[g], r["mutation"])] = float(r["ddG_pred"])
    out = []
    for r in csv.DictReader(open(var_csv, encoding="utf-8-sig")):
        g, v, pos = r["gene"], r["variante"], int(r["posicao"])
        if r["rotulo"] != "U" or not (DOM[g][0] <= pos <= DOM[g][1]):
            continue
        d, dd = dist[(g, pos)], ddg.get((g, pos, v[-1]))
        regra = d < 4 or (dd is not None and dd >= corte)
        sens, fp, lrp, lrn = calib[g]
        out.append({"gene": g, "variante": v, "clinvar": r["clinvar"], "gnomad_ac": r["gnomad_ac"],
                    "dist_dna": round(d, 1), "ddg_thermompnn": round(dd, 2) if dd is not None else "",
                    "alphamissense": r["alphamissense"], "regra": "positiva" if regra else "negativa",
                    "LR": round(lrp if regra else lrn, 2)})
    with open(saida, "w", newline="", encoding="utf-8-sig") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)
    for g in INI:
        s = [o for o in out if o["gene"] == g]
        pos_ = [o for o in s if o["regra"] == "positiva"]
        print(f"\n{g}: {len(s)} incertas no domínio | regra positiva {len(pos_)} | negativa {len(s) - len(pos_)}")
        no_clinvar = [o for o in s if o["clinvar"]]
        print(f"   das {len(no_clinvar)} incertas com registro no ClinVar:")
        for o in sorted(no_clinvar, key=lambda o: (o["regra"] != "positiva", -float(o["alphamissense"]))):
            print(f"   {o['variante']:7s} {o['regra']:9s} LR={o['LR']:<5} dist={o['dist_dna']:5} ΔΔG={o['ddg_thermompnn']:<5} "
                  f"AM={float(o['alphamissense']):.2f} gnomAD_AC={o['gnomad_ac']} | {o['clinvar'][:40]}")


if __name__ == "__main__":
    main()
