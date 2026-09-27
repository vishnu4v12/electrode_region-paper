"""Parse `mutool trace` XML into stroked sub-paths and text glyphs, in page coords (y down)."""
import xml.etree.ElementTree as ET
import numpy as np


def _tf(t):
    a, b, c, d, e, f = map(float, t.split())
    return lambda x, y: (a * x + c * y + e, b * x + d * y + f)


def parse(path):
    root = ET.parse(path).getroot()
    page = root.find("page")
    H = float(page.get("mediabox").split()[3])
    paths, glyphs = [], []
    for el in page.iter():
        if el.tag == "stroke_path":
            T = _tf(el.get("transform"))
            col = tuple(float(c) for c in el.get("color").split())
            lw = float(el.get("linewidth"))
            cur = []
            for cmd in el:
                if cmd.tag == "moveto":
                    if len(cur) > 1:
                        paths.append(dict(color=col, lw=lw, xy=np.array(cur)))
                    cur = [T(float(cmd.get("x")), float(cmd.get("y")))]
                elif cmd.tag == "lineto":
                    cur.append(T(float(cmd.get("x")), float(cmd.get("y"))))
                elif cmd.tag == "curveto":
                    cur.append(T(float(cmd.get("x3")), float(cmd.get("y3"))))
            if len(cur) > 1:
                paths.append(dict(color=col, lw=lw, xy=np.array(cur)))
        elif el.tag == "fill_text":
            T = _tf(el.get("transform"))
            for span in el.findall("span"):
                size = float(span.get("trm").split()[0])
                for g in span.findall("g"):
                    x, y = T(float(g.get("x")), float(g.get("y")))
                    glyphs.append((g.get("unicode"), x, y, size))
    # transform already has the -1 flip: y is measured from the top.
    return paths, glyphs, H
