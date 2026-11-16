"""Gera as figuras 1–5 do manuscrito (PNG 300 dpi + PDF) em agenesia/figuras/.

Uso: <venv-python> figuras.py <pasta_agenesia>
Regra do artigo: contato (<4 Å) OU ΔΔG >= 1,5 kcal/mol (corte da validação cruzada).
Cores (validadas para daltonismo): contato = azul, desestabiliza = laranja, nenhum/controle = cinza, patogênica = violeta.
"""
import csv
import math
import os
import re
import runpy
import sys
import zipfile
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

AZUL, LARANJA, CINZA, VIOLETA = "#2a78d6", "#eb6834", "#8a8984", "#6250d6"
TINTA, TINTA2, GRADE = "#0b0b0b", "#52514e", "#e6e5e1"
CORTE_DDG, CORTE_D = 1.5, 4.0
INI = {"PAX9": 6, "MSX1": 166}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "axes.edgecolor": TINTA2, "axes.labelcolor": TINTA, "xtick.color": TINTA2, "ytick.color": TINTA2,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "legend.frameon": False,
    "figure.facecolor": "white", "axes.facecolor": "white", "savefig.dpi": 300,
})


def ler(p):
    return list(csv.DictReader(open(p, encoding="utf-8-sig")))


def mecanismo(d, g):
    if d < CORTE_D:
        return "contato"
    if g is not None and g >= CORTE_DDG:
        return "desestab"
    return "nenhum"


COR = {"contato": AZUL, "desestab": LARANJA, "nenhum": CINZA}
NOME = {"contato": "DNA contact (<4 Å)", "desestab": f"Destabilizing (ΔΔG ≥ {CORTE_DDG})", "nenhum": "Neither"}


def salvar(fig, pasta, nome):
    fig.savefig(os.path.join(pasta, nome + ".png"), bbox_inches="tight")
    fig.savefig(os.path.join(pasta, nome + ".pdf"), bbox_inches="tight")
    plt.close(fig)


def fig1(base, out, res):
    """Distância ao DNA ao longo do domínio, com os sítios patogênicos coloridos pelo mecanismo."""
    dm = runpy.run_path(os.path.join(base, "scripts", "distancia_dna.py"))
    z = zipfile.ZipFile(os.path.join(base, "dados_brutos", "afserver", "lote_parcial.zip"))
    fig, axs = plt.subplots(2, 1, figsize=(7.2, 4.4), gridspec_kw={"height_ratios": [1.6, 1]})
    for ax, g in zip(axs, ("PAX9", "MSX1")):
        prot, dna = dm["ler_modelo"](z, dm["WT_DIR"][g])
        xs, ds = [], []
        for i in sorted(prot):
            xs.append(i + INI[g] - 1)
            ds.append(min(math.dist(a, b) for a, _ in prot[i] for b in dna))
        pts = [r for r in res if r["gene"] == g and r["variante"] != "WT"]
        vistos = {}
        for r in pts:
            pos = int(r["variante"][1:-1])
            mec = mecanismo(float(r["dist_dna_wt"]), float(r["ddg_thermompnn"]))
            vistos.setdefault(pos, []).append(mec)
        # barras: cinza claro = resíduo sem variante patogênica; sítios patogênicos em cinza escuro
        cores_barra = [TINTA2 if x in vistos else GRADE for x in xs]
        ax.bar(xs, ds, width=0.8, color=cores_barra, zorder=1)
        ax.axhline(CORTE_D, color=TINTA2, lw=0.6, ls="--", zorder=2)
        ax.text(xs[-1] + 0.5, CORTE_D, " 4 Å", va="center", ha="left", color=TINTA2, fontsize=7)
        # faixa superior: uma marca por variante, colorida pelo mecanismo (empilhadas no mesmo sítio)
        topo = max(ds) + 1.5
        for pos, mecs in vistos.items():
            for k, mec in enumerate(mecs):
                ax.scatter(pos, topo + 1.4 * k, s=20, color=COR[mec], edgecolor="white", linewidth=0.7, zorder=3)
        ax.set_ylim(0, topo + 1.4 * max(len(m) for m in vistos.values()) + 0.5)
        ax.set_xlim(xs[0] - 1, xs[-1] + 4)
        ax.set_ylabel("Distance to DNA (Å)")
        dom = "paired domain" if g == "PAX9" else "homeodomain"
        ax.set_title(f"{'A' if g == 'PAX9' else 'B'}   {g} {dom}: distance of each residue to DNA (AlphaFold 3)", loc="left")
        ax.yaxis.grid(True, color=GRADE, lw=0.5)
        ax.set_axisbelow(True)
    axs[1].set_xlabel("Residue (UniProt numbering)")
    for mec in ("contato", "desestab", "nenhum"):
        axs[0].scatter([], [], s=22, color=COR[mec], label=NOME[mec])
    from matplotlib.patches import Patch
    hs, ls = axs[0].get_legend_handles_labels()
    hs.append(Patch(color=TINTA2))
    ls.append("Residue with pathogenic variant (bar)")
    axs[0].legend(hs, ls, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=4, fontsize=7, handletextpad=0.3)
    fig.tight_layout()
    salvar(fig, out, "Fig1_distancia_dna")


def fig2(out, res):
    """Distância ao DNA x ΔΔG: os dois mecanismos."""
    fig, ax = plt.subplots(figsize=(4.6, 3.8))
    from matplotlib.patches import Rectangle
    ax.add_patch(Rectangle((0, -0.6), CORTE_D, 4.2, color=AZUL, alpha=0.06, lw=0))
    ax.add_patch(Rectangle((CORTE_D, CORTE_DDG), 16.5 - CORTE_D, 3.6 - CORTE_DDG, color=LARANJA, alpha=0.06, lw=0))
    ax.axvline(CORTE_D, color=TINTA2, lw=0.6, ls="--")
    ax.axhline(CORTE_DDG, color=TINTA2, lw=0.6, ls="--")
    # rótulos com deslocamento manual (pontos), para não sobrepor
    rot = {("PAX9", "R28P"): (5, 2), ("PAX9", "A108P"): (-8, 6), ("PAX9", "L27P"): (5, 2), ("PAX9", "I56N"): (4, 5),
           ("PAX9", "R26W"): (6, 2), ("PAX9", "G51S"): (6, -3), ("PAX9", "L21P"): (6, 4), ("PAX9", "G6R"): (5, 3),
           ("PAX9", "R47W"): (5, 3), ("MSX1", "L211R"): (6, -9), ("MSX1", "R202P"): (-5, 6), ("MSX1", "R202H"): (-6, -2)}
    for r in res:
        if r["variante"] == "WT":
            continue
        d, g = float(r["dist_dna_wt"]), float(r["ddg_thermompnn"])
        mec = mecanismo(d, g)
        marker = "o" if r["gene"] == "PAX9" else "s"
        ax.scatter(d, g, s=26, marker=marker, color=COR[mec], edgecolor="white", linewidth=0.8, zorder=3)
        off = rot.get((r["gene"], r["variante"]))
        if off:
            lab = r["variante"] if r["gene"] == "PAX9" else f"MSX1 {r['variante']}"
            ax.annotate(lab, (d, g), xytext=off, textcoords="offset points", fontsize=6.5, color=TINTA,
                        ha="left" if off[0] >= 0 else "right")
    ax.set_xlim(0, 16.5)
    ax.set_ylim(-0.6, 3.6)
    ax.set_xlabel("Distance of mutated residue to DNA (Å)")
    ax.set_ylabel("Predicted ΔΔG folding (kcal/mol, ThermoMPNN)")
    ax.text(0.2, 3.45, "DNA contact", color=AZUL, fontsize=7, va="top")
    ax.text(16.2, 3.45, "Destabilizing", color=LARANJA, fontsize=7, va="top", ha="right")
    for mec in ("contato", "desestab", "nenhum"):
        ax.scatter([], [], s=22, color=COR[mec], label=NOME[mec])
    ax.scatter([], [], s=22, marker="o", color=TINTA2, label="PAX9")
    ax.scatter([], [], s=22, marker="s", color=TINTA2, label="MSX1")
    ax.legend(loc="lower right", fontsize=6.5, handletextpad=0.3)
    ax.set_title("Pathogenic missense variants: two mechanisms", loc="left")
    salvar(fig, out, "Fig2_mecanismos")


def fig3(out, ctrl):
    """Patogênicas x controle gnomAD (PAX9): ΔΔG e proporções."""
    P = [r for r in ctrl if r["gene"] == "PAX9" and r["grupo"] == "patogenica"]
    C = [r for r in ctrl if r["gene"] == "PAX9" and r["grupo"] == "controle_gnomad" and int(r["gnomad_ac"]) >= 2]
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.0, 3.2), gridspec_kw={"width_ratios": [1, 1.3]})
    import random
    random.seed(1)
    for i, (grp, cor, nome) in enumerate(((P, VIOLETA, f"Pathogenic\n(n={len(P)})"), (C, CINZA, f"gnomAD AC≥2\n(n={len(C)})"))):
        ys = [float(r["ddg_thermompnn"]) for r in grp if r["ddg_thermompnn"] != ""]
        xs = [i + random.uniform(-0.18, 0.18) for _ in ys]
        a.scatter(xs, ys, s=12, color=cor, edgecolor="white", linewidth=0.5, zorder=3)
        m = median(ys)
        a.plot([i - 0.28, i + 0.28], [m, m], color=TINTA, lw=1.6, zorder=4, solid_capstyle="round")
        a.text(i + 0.31, m, f"{m:.2f}", va="center", fontsize=7, color=TINTA)
    a.axhline(CORTE_DDG, color=TINTA2, lw=0.6, ls="--")
    a.set_xticks([0, 1], [f"Pathogenic\n(n={len(P)})", f"gnomAD AC≥2\n(n={len(C)})"])
    a.set_xlim(-0.6, 1.7)
    a.set_ylabel("Predicted ΔΔG (kcal/mol)")
    a.set_title("A   PAX9 domain stability (p = 4×10⁻⁵)", loc="left")
    a.yaxis.grid(True, color=GRADE, lw=0.5)
    a.set_axisbelow(True)

    def frac(grp, f):
        return 100 * sum(map(f, grp)) / len(grp)
    contato = lambda r: float(r["dist_dna"]) < CORTE_D
    desest = lambda r: r["ddg_thermompnn"] != "" and float(r["ddg_thermompnn"]) >= CORTE_DDG
    regra = lambda r: contato(r) or desest(r)
    cats = ["DNA contact", f"ΔΔG ≥ {CORTE_DDG}", "Contact or\nΔΔG ≥ 1.5"]
    vp = [frac(P, contato), frac(P, desest), frac(P, regra)]
    vc = [frac(C, contato), frac(C, desest), frac(C, regra)]
    w = 0.36
    for k in range(3):
        b.bar(k - w / 2 - 0.01, vp[k], width=w, color=VIOLETA, zorder=3)
        b.bar(k + w / 2 + 0.01, vc[k], width=w, color=CINZA, zorder=3)
        b.text(k - w / 2 - 0.01, vp[k] + 1.5, f"{vp[k]:.0f}%", ha="center", fontsize=7, color=TINTA)
        b.text(k + w / 2 + 0.01, vc[k] + 1.5, f"{vc[k]:.0f}%", ha="center", fontsize=7, color=TINTA)
    b.set_xticks(range(3), cats)
    b.set_ylim(0, 100)
    b.set_ylabel("Variants (%)")
    b.set_title("B   Pathogenic vs population variants (PAX9)", loc="left")
    from matplotlib.patches import Patch
    b.legend(handles=[Patch(color=VIOLETA, label="Pathogenic"), Patch(color=CINZA, label="gnomAD AC≥2")],
             loc="upper left", fontsize=7)
    b.yaxis.grid(True, color=GRADE, lw=0.5)
    b.set_axisbelow(True)
    fig.tight_layout()
    salvar(fig, out, "Fig3_controle_gnomad")


def fig4(base, out, res):
    """Dentes ausentes por mecanismo e variantes nulas."""
    AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E", "Gly": "G",
           "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S",
           "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
    mis, nul = {}, []
    for ln in open(os.path.join(base, "dados_brutos", "literatura", "intarak2023_supl.txt"), encoding="utf-8"):
        c = [x.strip() for x in ln.split("|")]
        if len(c) < 6:
            continue
        try:
            dentes = float(c[4])
        except ValueError:
            continue
        m = re.fullmatch(r"([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})", c[2])
        if m and m.group(1) in AA3 and m.group(3) in AA3:
            mis.setdefault(AA3[m.group(1)] + m.group(2) + AA3[m.group(3)], dentes)
        elif "fs" in c[2] or c[2].endswith("*"):
            nul.append(dentes)
    rmap = {r["variante"]: r for r in res if r["gene"] == "PAX9"}
    grupos = {"contato": [], "desestab": [], "nenhum": []}
    for v, d in mis.items():
        r = rmap.get(v)
        if r:
            grupos[mecanismo(float(r["dist_dna_wt"]), float(r["ddg_thermompnn"]))].append(d)
    import random
    random.seed(2)
    fig, ax = plt.subplots(figsize=(4.6, 3.4))
    series = [("DNA contact", grupos["contato"], AZUL), (f"Destabilizing", grupos["desestab"], LARANJA),
              ("Neither", grupos["nenhum"], CINZA), ("Null variants\n(frameshift/nonsense)", nul, TINTA2)]
    for i, (nome, ys, cor) in enumerate(series):
        if not ys:
            continue
        xs = [i + random.uniform(-0.17, 0.17) for _ in ys]
        ax.scatter(xs, ys, s=14, color=cor, edgecolor="white", linewidth=0.5, zorder=3)
        m = median(ys)
        ax.plot([i - 0.27, i + 0.27], [m, m], color=TINTA, lw=1.6, zorder=4, solid_capstyle="round")
        ax.text(i + 0.3, m, f"{m:g}", va="center", fontsize=7, color=TINTA)
    ax.set_xticks(range(4), [f"{s[0]}\n(n={len(s[1])})" for s in series], fontsize=7)
    ax.set_ylabel("Missing teeth (excl. third molars)")
    ax.set_title("PAX9: severity by mechanism (missense vs null p = 0.001)", loc="left")
    ax.yaxis.grid(True, color=GRADE, lw=0.5)
    ax.set_axisbelow(True)
    salvar(fig, out, "Fig4_gravidade")


def fig5(out, comp):
    """AlphaGenome 2 x SpliceAI para variantes não codificantes."""
    fig, ax = plt.subplots(figsize=(4.6, 3.8))
    cores = {"P": VIOLETA, "B": CINZA, "U": TINTA}
    nomes = {"P": "Pathogenic", "B": "Benign", "U": "VUS"}
    # rótulo único para os dois VUS do MSX1 (mesmo ponto) e para o grupo de VUS do PAX9
    # posição do texto em fração dos eixos (área vazia), com linha até o ponto
    destaque = {"c.470-9G>A": ("MSX1 c.470-9G>A and\nc.470-10C>A (VUS)", (0.62, 0.80)),
                "c.771+3G>T": ("PAX9 splice-site VUS\n(c.771+1/+3/+5, c.772-2/-3)", (0.62, 0.62))}
    for g in ("B", "U", "P"):
        sub = [r for r in comp if r["grupo"] == g and r["spliceai_max"] not in ("", None)]
        xs = [max(float(r["spliceai_max"]), 0.005) for r in sub]
        ys = [max(float(r["max_splice_sites_q"]), 0.001) for r in sub]
        if g == "U":
            ax.scatter(xs, ys, s=20, facecolor="white", edgecolor=TINTA, linewidth=0.8, zorder=3, label=nomes[g])
        else:
            ax.scatter(xs, ys, s=20, color=cores[g], edgecolor="white", linewidth=0.6, zorder=3, label=nomes[g])
        for r, x, y in zip(sub, xs, ys):
            for k, (lab, off) in destaque.items():
                if r["title"].endswith(k):
                    ax.annotate(lab, (x, y), xytext=off, textcoords="axes fraction", ha="right", va="center",
                                fontsize=6.3, color=TINTA,
                                arrowprops=dict(arrowstyle="-", color=TINTA2, lw=0.5, shrinkA=0, shrinkB=3))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("SpliceAI (max Δ score; values < 0.005 shown at 0.005)")
    ax.set_ylabel("AlphaGenome 2 splice-site score (max, 3 cell types)")
    ax.set_title("Non-coding variants: AlphaGenome vs SpliceAI (AUC 1.0 both)", loc="left")
    ax.legend(loc="upper left", fontsize=7)
    ax.grid(True, color=GRADE, lw=0.5, which="major")
    ax.set_axisbelow(True)
    salvar(fig, out, "Fig5_alphagenome_spliceai")


def main():
    base = sys.argv[1]
    out = os.path.join(base, "figuras")
    os.makedirs(out, exist_ok=True)
    res = [r for r in ler(os.path.join(base, "resultados_afserver.csv"))]
    ctrl = ler(os.path.join(base, "grupo_controle.csv"))
    comp = ler(os.path.join(base, "comparacao_ag1_ag2_spliceai.csv"))
    fig1(base, out, res)
    fig2(out, res)
    fig3(out, ctrl)
    fig4(base, out, res)
    fig5(out, comp)
    print("figuras em", out, sorted(os.listdir(out)))


if __name__ == "__main__":
    main()
