"""
Build the Editorial Manager upload set in ../submission/: a single self-contained
manuscript.tex (all \\input files and the compiled bibliography inlined, figures
renamed Figure_N.pdf in order of appearance), the figure files, the highlights as
a .docx, and the compiled manuscript and cover letter PDFs.
Run after compiling main.tex (main.bbl must be current).
"""
import os
import re
import shutil

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "submission")


def read(p):
    return open(os.path.join(HERE, p), encoding="utf-8").read()


def inline(text):
    def rep(m):
        name = m.group(1)
        path = name if name.endswith(".tex") else name + ".tex"
        return inline(read(path).rstrip("\n"))
    return re.sub(r"\\input\{([^}]+)\}", rep, text)


def main():
    os.makedirs(OUT, exist_ok=True)
    tex = inline(read("main.tex"))
    tex = tex.replace("\\bibliographystyle{elsarticle-harv}\n", "")
    tex = re.sub(r"\\bibliography\{references\}", lambda m: read("main.bbl").rstrip("\n"), tex)
    figs = []
    for f in re.findall(r"\\includegraphics(?:\[[^]]*\])?\{figures/([^}]+)\}", tex):
        if f not in figs:
            figs.append(f)
    for i, f in enumerate(figs, 1):
        tex = tex.replace("{figures/" + f + "}", "{Figure_%d.pdf}" % i)
        shutil.copy(os.path.join(HERE, "figures", f), os.path.join(OUT, "Figure_%d.pdf" % i))
    open(os.path.join(OUT, "manuscript.tex"), "w", encoding="utf-8").write(tex)

    import docx
    d = docx.Document()
    d.add_heading("Highlights", level=1)
    for line in read("sections/highlights.tex").splitlines():
        s = line.replace("\\item ", "").strip()
        if s:
            assert len(s) <= 85, s
            d.add_paragraph(s, style="List Bullet")
    d.save(os.path.join(OUT, "highlights.docx"))
    shutil.copy(os.path.join(HERE, "main.pdf"), os.path.join(OUT, "manuscript_compiled.pdf"))
    shutil.copy(os.path.join(HERE, "cover_letter.pdf"), os.path.join(OUT, "cover_letter.pdf"))
    print("figures:", figs)


if __name__ == "__main__":
    main()
