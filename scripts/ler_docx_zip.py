"""Extrai o texto das tabelas de arquivos .docx dentro de um .zip (sem executar nada do arquivo).

Uso: python -I ler_docx_zip.py <arquivo.zip> <saida.txt>
"""
import io
import re
import sys
import zipfile
from xml.etree import ElementTree as ET

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"


def docx_text(data):
    with zipfile.ZipFile(io.BytesIO(data)) as z:
        root = ET.fromstring(z.read("word/document.xml"))
    out = []
    for el in root.iter():
        if el.tag == W + "tbl":
            for tr in el.iter(W + "tr"):
                cells = ["".join(t.text or "" for t in tc.iter(W + "t")) for tc in tr.iter(W + "tc")]
                out.append(" | ".join(c.strip() for c in cells))
            out.append("")
    return "\n".join(out)


def main():
    zpath, saida = sys.argv[1], sys.argv[2]
    partes = []
    with zipfile.ZipFile(zpath) as z:
        for n in z.namelist():
            if n.lower().endswith(".docx"):
                partes.append(f"### {n}\n" + docx_text(z.read(n)))
            else:
                partes.append(f"### {n} (ignorado)")
    texto = "\n".join(partes)
    with open(saida, "w", encoding="utf-8") as f:
        f.write(texto)
    print(texto[:200])
    print("missense p.:", sorted(set(re.findall(r"p\.\(?([A-Z][a-z]{2}\d+[A-Z][a-z]{2})\)?", texto))))


if __name__ == "__main__":
    main()
