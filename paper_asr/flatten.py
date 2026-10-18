"""
Build the Editorial Manager upload set in ../submission/ (run after compiling main.tex):

  source/manuscript.tex     single self-contained LaTeX file (all \\input files and the
                            compiled bibliography inlined; figures named Figure_N.pdf)
  source/Figure_N.pdf       figures in order of appearance (also uploaded as figure items)
  source/references.bib     BibTeX database (for the typesetter)
  source/elsarticle.cls, source/elsarticle-harv.bst   Elsevier class and style
  highlights.docx, highlights.txt, abstract.txt, cover_letter.pdf, manuscript_compiled.pdf
"""
import os
import re
import shutil
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "submission")
SRC = os.path.join(OUT, "source")


def read(p):
    return open(os.path.join(HERE, p), encoding="utf-8").read()


def inline(text):
    def rep(m):
        name = m.group(1)
        path = name if name.endswith(".tex") else name + ".tex"
        return inline(read(path).rstrip("\n"))
    return re.sub(r"\\input\{([^}]+)\}", rep, text)


def plain(tex, macros):
    """LaTeX abstract -> plain UTF-8 text for the submission form."""
    for k in sorted(macros, key=len, reverse=True):
        tex = tex.replace("\\" + k + "{}", macros[k]).replace("\\" + k, macros[k])
    rep = [(r"\SI{8}{s}", "8 s"), (r"\textdegree{}", "°"), (r"\textdegree", "°"), ("K$\\alpha$", "Kα"),
           ("$-$", "−"), ("{,}", ","), ("--", "–"), ("\\%", "%")]
    for a, b in rep:
        tex = tex.replace(a, b)
    tex = re.sub(r"\$r=([^$]+)\$", r"r = \1", tex)
    tex = tex.replace("$", "")
    assert "\\" not in tex, tex[tex.index("\\") - 30:tex.index("\\") + 30]
    return re.sub(r"\s+", " ", tex).strip()


def kpse(name):
    return subprocess.run(["kpsewhich", name], capture_output=True, text=True).stdout.strip()


def main():
    os.makedirs(SRC, exist_ok=True)
    tex = inline(read("main.tex"))
    tex = tex.replace("\\bibliographystyle{elsarticle-harv}\n", "")
    tex = re.sub(r"\\bibliography\{references\}", lambda m: read("main.bbl").rstrip("\n"), tex)
    figs = []
    for f in re.findall(r"\\includegraphics(?:\[[^]]*\])?\{figures/([^}]+)\}", tex):
        if f not in figs:
            figs.append(f)
    for i, f in enumerate(figs, 1):
        tex = tex.replace("{figures/" + f + "}", "{Figure_%d.pdf}" % i)
        shutil.copy(os.path.join(HERE, "figures", f), os.path.join(SRC, "Figure_%d.pdf" % i))
    open(os.path.join(SRC, "manuscript.tex"), "w", encoding="utf-8").write(tex)
    shutil.copy(os.path.join(HERE, "references.bib"), os.path.join(SRC, "references.bib"))
    # bibliography is already inlined in manuscript.tex; the .bbl is supplied as well for the typesetter
    shutil.copy(os.path.join(HERE, "main.bbl"), os.path.join(SRC, "manuscript.bbl"))
    for f in ["elsarticle.cls", "elsarticle-harv.bst"]:
        p = kpse(f)
        if p:
            shutil.copy(p, os.path.join(SRC, f))

    hl = [ln.replace("\\item ", "").strip() for ln in read("sections/highlights.tex").splitlines() if ln.strip()]
    for s in hl:
        assert len(s) <= 85, (len(s), s)
    open(os.path.join(OUT, "highlights.txt"), "w", encoding="utf-8").write("\n".join("• " + s for s in hl) + "\n")
    import docx
    d = docx.Document()
    d.add_heading("Highlights", level=1)
    for s in hl:
        d.add_paragraph(s, style="List Bullet")
    d.save(os.path.join(OUT, "highlights.docx"))

    macros = dict(re.findall(r"\\newcommand\{\\(\w+)\}\{(.*)\}", read("numbers.tex")))
    ab = plain(read("sections/abstract.tex"), macros)
    open(os.path.join(OUT, "abstract.txt"), "w", encoding="utf-8").write(ab + "\n")

    shutil.copy(os.path.join(HERE, "main.pdf"), os.path.join(OUT, "manuscript_compiled.pdf"))
    shutil.copy(os.path.join(HERE, "cover_letter.pdf"), os.path.join(OUT, "cover_letter.pdf"))
    print("figures:", figs)
    print("highlights (chars):", [len(s) for s in hl])
    print("abstract words:", len(ab.split()))


if __name__ == "__main__":
    main()
