"""Render two hierarchy/boundary variants of the problem-framing figure (based on the editorial layout)."""
import sys
import os
from pathlib import Path

from Bio.PDB import PDBParser
import fitz
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import (Circle, FancyArrowPatch, FancyBboxPatch,
                                PathPatch, Rectangle)
from matplotlib.path import Path as DrawPath
import numpy as np
from PIL import Image
from scipy.interpolate import CubicSpline

HERE = Path(__file__).resolve().parent
ROOT = HERE
OUT = HERE / "rendered"
MODE = "glossy"
FLAT = os.environ.get("FLAT") == "1"
GRAYSUP = os.environ.get("GRAYSUP") == "1"
PASTEL = os.environ.get("PASTEL") == "1"

MAIN = MODE == "mainsec"
SUFFIX = "flat" if FLAT else MODE
if FLAT and GRAYSUP:
    SUFFIX = "flat_hero"
if FLAT and PASTEL:
    SUFFIX = "flat_pastel"
DEST = OUT / "variants" / SUFFIX
DEST.mkdir(parents=True, exist_ok=True)

parser = PDBParser(QUIET=True)
PDB_DIR = ROOT / "outputs/length_research_handoff_20260926/DS04_work/sources/PRA310_official/PRA310/PDBs"
SOURCES = {
    "2C4R": PDB_DIR / "2C4R.pdb",                     # RNase E-ssRNA (landscape L01)
    "4QVC": PDB_DIR / "4QVC.pdb",                     # Hfq-OxyS (landscape L10/L20)
    "2ZZN": PDB_DIR / "2ZZN.pdb",                     # Trm5-tRNA Cys (landscape L13)
    "4NGD": ROOT / "outputs/dicer_input_preparation_20260923/inputs/structures/4NGD.pdb",  # Dicer-siRNA
}
CHAINS = {"2C4R": ("L", "R"), "4QVC": ("ABCDEF", "G"), "2ZZN": ("A", "C"), "4NGD": ("A", "B")}
W, H = 284, 141
INK, GRAY, RULE = "#172735", "#536270", "#CCD5DC"
BLUE, ORANGE = "#00527F", "#A74900"
DARK = "#1C4048"
TEAL = "#427C7B"
CARD_FC, CARD_EC = "#F8FBFC", "#D2DFE3"
BAND_FC, BAND_LINE = "#EDF4F7", "#C6D6DC"
mpl.rcParams.update({"font.family": "Arial", "mathtext.fontset": "dejavusans",
                     "pdf.fonttype": 42, "svg.fonttype": "none"})
fig = plt.figure(figsize=(W / 25.4, H / 25.4), facecolor="white")
ax = fig.add_axes([0, 0, 1, 1], xlim=(0, W), ylim=(H, 0))
ax.set_axis_off()


def txt(x, y, s, *, c=INK, size=11.5, ha="left", bold=False):
    ax.text(x, y, s, va="center", ha=ha, color=c, fontsize=size,
            fontweight="bold" if bold else "normal", zorder=15)


def line(x, y, color=RULE, lw=.8, style="-", z=2):
    ax.plot(x, y, color=color, lw=lw, ls=style, zorder=z, solid_capstyle="round")


def dot(x, y, marker="o", color=BLUE, size=5.5, fill=True, z=10):
    ax.plot(x, y, ls="none", marker=marker, ms=size, mec=color,
            mfc=color if fill else "white", mew=1.15, zorder=z)


def model(pdb):
    return parser.get_structure(pdb, SOURCES[pdb])[0]


def backbone(pdb, m):
    paths, bases = [], []
    protein, rna = CHAINS[pdb]
    for chain in protein + rna:
        if chain not in m:
            continue
        is_rna = chain in rna
        names = ("C4'", "C4*") if is_rna else ("CA",)
        residues = [r for r in m[chain] if r.id[0] == " "]
        trace = []
        for r in residues:
            atom = next((r[name] for name in names if name in r), None)
            if atom is None:
                continue
            trace.append(atom.coord)
            if is_rna:
                base = next((r[name] for name in ("N9", "N1") if name in r), None)
                if base is not None:
                    bases.append((atom.coord, base.coord))
        if len(trace) > 1:
            paths.append((np.asarray(trace), is_rna))
    assert any(is_rna for _, is_rna in paths)
    assert any(not is_rna for _, is_rna in paths)
    return paths, bases


def structure(pdb, x, y, span, active, rainbow_style=True):
    """RFdiffusion-inspired rendering: rainbow N->C traces with white halo."""
    m = model(pdb)
    paths, bases = backbone(pdb, m)
    cloud = np.concatenate([trace for trace, _ in paths])
    center = cloud.mean(axis=0)
    axes = np.linalg.svd(cloud - center, full_matrices=False)[2]
    projected = (cloud - center) @ axes.T
    extent = np.ptp(projected[:, :2], axis=0).max()
    offset = (projected.max(axis=0) + projected.min(axis=0)) / 2

    def project(points):
        return (((np.asarray(points) - center) @ axes.T) - offset) * span / extent + [x, y, 0]

    segments = []
    for trace, is_rna in paths:
        n = len(trace)
        gaps = np.flatnonzero(np.linalg.norm(np.diff(trace, axis=0), axis=1) >
                              (9 if is_rna else 5)) + 1
        bounds = np.concatenate(([0], gaps, [n]))
        for a, b in zip(bounds[:-1], bounds[1:]):
            part = trace[a:b]
            if len(part) < 3:
                continue
            sm = CubicSpline(np.arange(len(part)), project(part))(
                np.linspace(0, len(part) - 1, (len(part) - 1) * 5 + 1))
            fracs = np.linspace(a, b - 1, len(sm)) / max(n - 1, 1)
            for i, (pa, pb) in enumerate(zip(sm[:-1], sm[1:])):
                segments.append((.5 * (pa[2] + pb[2]), pa[:2], pb[:2], is_rna, float(fracs[i])))
    segments.sort(key=lambda item: item[0])
    lo, hi = segments[0][0], segments[-1][0]

    rainbow = mpl.colors.LinearSegmentedColormap.from_list(
        "rf", ["#3B5BA5", "#3F9BC0", "#57A96B", "#DFC144", "#D97B3C", "#B84A3F"])
    colors, widths, coords = [], [], []
    for depth, pa, pb, is_rna, frac in segments:
        strength = (depth - lo) / (hi - lo) if hi > lo else 0.5
        shade = 0.80 + 0.30 * strength
        if is_rna:
            rgb = np.array(mpl.colors.to_rgb("#00527F"))
        elif rainbow_style:
            rgb = np.array(rainbow(np.clip(frac, 0, 1)))
            if PASTEL and not active:
                rgb = 0.62 * rgb + 0.38
        else:
            rgb = np.array(mpl.colors.to_rgb("#9AAEB8"))
        colors.append(np.clip(rgb * shade, 0, 1))
        widths.append((2.8 if is_rna else 2.2) if active else (2.4 if is_rna else 1.85))
        coords.append([pa, pb])
    ax.add_collection(LineCollection(coords, colors="#FFFFFF",
                                     linewidths=[w + 1.1 for w in widths],
                                     capstyle="round", joinstyle="round", zorder=4.5))
    ax.add_collection(LineCollection(coords, colors=colors, linewidths=widths,
                                     capstyle="round", joinstyle="round", zorder=5))
    for atoms in bases:
        a, b = project(atoms)
        line([a[0], b[0]], [a[1], b[1]], "#124E6E", .7 if active else .55, z=6)
    return m, project


def interface_pair(m):
    protein = [r for r in m["B"] if r.id[0] == " " and "CA" in r]
    rna = [r for r in m["Q"] if r.id[0] == " "]
    candidates = []
    for pr in protein:
        side = [a for a in pr if a.name not in {"N", "CA", "C", "O"} and a.element != "H"]
        for nu in rna:
            for a in side:
                for b in nu:
                    if a.element in {"N", "O"} and b.element in {"N", "O"}:
                        distance = np.linalg.norm(a.coord - b.coord)
                        if 2.5 <= distance <= 3.5:
                            candidates.append((distance, pr, nu, a, b))
    assert candidates
    return min(candidates, key=lambda item: item[0])


def draw_contact(pair, x, y, span):
    _, protein, rna, contact_a, contact_b = pair
    atoms = [[a for a in r if a.element != "H"] for r in (protein, rna)]
    cloud = np.array([a.coord for group in atoms for a in group])
    center = cloud.mean(axis=0)
    axes = np.linalg.svd(cloud - center, full_matrices=False)[2]
    cloud = (cloud - center) @ axes[:2].T
    scale = span / np.ptp(cloud, axis=0).max()
    positions = {id(a): (a.coord - center) @ axes[:2].T * scale + [x, y]
                 for group in atoms for a in group}
    a, b = positions[id(contact_a)], positions[id(contact_b)]
    line([a[0], b[0]], [a[1], b[1]], GRAY, .9, "--", 7)
    for group, color in zip(atoms, ("#6A7883", BLUE)):
        for i, a in enumerate(group):
            for b in group[i + 1:]:
                threshold = 1.95 if "P" in (a.element, b.element) else 1.78
                if .9 < np.linalg.norm(a.coord - b.coord) < threshold:
                    pa, pb = positions[id(a)], positions[id(b)]
                    line([pa[0], pb[0]], [pa[1], pb[1]], color, 1.45, z=8)
            px, py = positions[id(a)]
            dot(px, py, color=color, size=2.7 if a.element == "C" else 3.4,
                fill=a.element == "C", z=9)


def curve(points, c=TEAL, lw=1, z=4):
    ax.add_patch(PathPatch(DrawPath(points, [DrawPath.MOVETO] + [DrawPath.CURVE4] * (len(points) - 1)),
                           facecolor="none", edgecolor=c, lw=lw, capstyle="round", zorder=z))


def arrow(x0, y0, x1, y1, c="#8FAEB2", lw=.9):
    ax.add_patch(FancyArrowPatch((x0, y0), (x1, y1), arrowstyle="-|>", lw=lw,
                                 mutation_scale=7.5, color=c, zorder=8,
                                 shrinkA=0, shrinkB=0))


def export(name):
    path = DEST / name
    for ext in ("png", "pdf", "svg"):
        fig.savefig(path.with_suffix(f".{ext}"), dpi=285, facecolor=fig.get_facecolor())
    with fitz.open(path.with_suffix(".pdf")) as doc:
        page = doc[0]
        for block in page.get_text("blocks"):
            assert page.rect.contains(fitz.Rect(block[:4])), block
        pix = page.get_pixmap(matrix=fitz.Matrix(2, 2), alpha=False)
        pix.save(DEST / f"{name}_preview.png")


RENDER = ROOT / "outputs" / "problem_framing_ribbon_20261107" / "renders"


CARD_CLIP = {}
CARD_CLIP["a"] = (2.0, 16.8, 60.6, 126.8)
CARD_CLIP["b1"] = (66.0, 16.8, 180.8, 68.9)
CARD_CLIP["c1"] = (188.0, 16.8, 280.2, 68.9)


def _clear_white(im, threshold=250):
    """White background -> transparent so the panel colour shows through (no visible image box)."""
    from PIL import Image as _Image
    arr = np.asarray(im.convert("RGB")).astype(np.uint8)
    alpha = np.where(np.min(arr, axis=2) >= threshold, 0, 255).astype(np.uint8)
    rgba = np.dstack([arr, alpha])
    return _Image.fromarray(rgba, "RGBA")



def _clip_patch(key):
    if key not in CARD_CLIP:
        return None
    x0, y0, x1, y1 = CARD_CLIP[key]
    patch = FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                           boxstyle="round,pad=0,rounding_size=2.0",
                           fc="none", ec="none", zorder=-10)
    ax.add_patch(patch)
    return patch


def struct_img(name, x, y, width, z=6, clip=None, bleed=0.0):
    """Place a PyMOL ribbon render (transparent background) with optional bleed + card clip."""
    from PIL import Image as _Image
    im = _Image.open(RENDER / name).convert("RGBA")
    arr = np.asarray(im)
    mask = arr[:, :, 3] > 8
    ys, xs = np.where(mask)
    if len(xs):
        im = im.crop((int(xs.min()), int(ys.min()), int(xs.max()) + 1, int(ys.max()) + 1))
        w0, h0 = im.size
        if w0 > 12 and h0 > 12:
            im = im.crop((4, 4, w0 - 4, h0 - 4))
    if bleed:
        from PIL import ImageOps as _Ops
        im = _Ops.expand(im, border=int(round(bleed * max(im.size))), fill=(0, 0, 0, 0))
    w, h = im.size
    full = float(width) * (1.0 + 2.0 * bleed)
    height = full * h / w
    ext = [x - full / 2, x + full / 2, y + height / 2, y - height / 2]
    art = ax.imshow(im, extent=ext, aspect="auto", zorder=z, interpolation="none")
    clip_patch = _clip_patch(clip)
    if clip_patch is not None:
        art.set_clip_path(clip_patch)


def ribbon_image():
    path = ROOT / "outputs/figure1_redesign_20261023/complex_1a9n.png"
    im = Image.open(path).convert("RGBA")
    return im.crop(im.getbbox())


def schematic(x, y, s=1, site=0):
    def pt(u, v):
        return x + s * u, y + s * v
    body = [(-6, -3), (-5, -7), (-1, -8), (1, -6), (5, -7), (7, -4), (6, -1),
            (5, 2), (2, 3), (-1, 2), (-4, 4), (-8, 1), (-6, -3)]
    path = DrawPath([pt(u, v) for u, v in body], [DrawPath.MOVETO] + [DrawPath.CURVE4] * 12)
    ax.add_patch(PathPatch(path, facecolor="#AABCC4", edgecolor="#647D85", lw=.85, zorder=4))
    strand = [(-7, 3), (-4, 5), (-1, 3), (1, 4), (3, 5), (4, 6), (7, 4)]
    curve([pt(u, v) for u, v in strand], BLUE, 1.4, 6)
    for u, v in [(-4, 4), (-1, 4), (2, 4.4), (5, 5.2)]:
        p, q = pt(u, v), pt(u + .4, v + 1.4)
        line([p[0], q[0]], [p[1], q[1]], BLUE, .65, z=5)
    for i, (u, v) in enumerate([(-4, -.1), (-1, -2), (2, -2.4), (4, .1)], 1):
        xx, yy = pt(u, v)
        dot(xx, yy, "D" if i == site else "o", ORANGE if i == site else "#56707B",
            4.6 if i == site else 2.25, fill=i == site)


def name(x, y, label, c=GRAY):
    txt(x, y, label, c=c, ha="center")


def cap(x, y, s):
    txt(x, y, s, c=DARK, bold=True, size=(10.4 if MAIN else 10.9) if MODE.startswith("nature") else (11.2 if MAIN else 11.8))


def level_marker(x, y):
    ax.add_patch(Rectangle((x - 3.1, y - 2.2), 1.15, 4.4, fc=BLUE, ec="none", zorder=6))


def token_row(x, y, letters):
    for i, letter in enumerate(letters):
        xx = x + i * 2.8
        ax.add_patch(Rectangle((xx, y - 1.35), 2.35, 2.7,
                              fc="#D9E7E8" if i != 2 else "#EBC6A9", ec="none", zorder=3))
        txt(xx + 1.18, y, letter, size=6.3, c=DARK, ha="center")


def graph(x, y):
    nodes = np.array([[-4, -2], [-1, -4], [3, -2], [4, 2], [0, 4], [-3, 1]])
    for a, b in [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0), (1, 5), (2, 4)]:
        line([x + nodes[a, 0], x + nodes[b, 0]],
             [y + nodes[a, 1], y + nodes[b, 1]], "#7C9CA2", .7, z=3)
    for i, (dx, dy) in enumerate(nodes):
        dot(x + dx, y + dy, "o", ORANGE if i == 2 else TEAL, 2.75)


# ---------------------------------------------------------------- page furniture
CARDS = {"a": (3.0, 15.8, 61.2, 126.4),
         "b": (65.4, 15.8, 181.0, 126.4),
         "c": (185.2, 15.8, 281.0, 126.4)}
PANELS = [(5, "a", "Molecular context"), (69, "b", "Two comparison scales"), (191, "c", "Research strategy")]

if MODE == "cards":
    for x0, y0, x1, y1 in CARDS.values():
        ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                    boxstyle="round,pad=0,rounding_size=1.6",
                                    fc=CARD_FC, ec=CARD_EC, lw=.8, zorder=.4))
    for x, tag, title in PANELS:
        ax.add_patch(FancyBboxPatch((x - 1.4, 2.6), 6.4, 6.4,
                                    boxstyle="round,pad=0,rounding_size=1.0",
                                    fc=BLUE, ec="none", zorder=6))
        txt(x + 1.8, 5.85, tag, c="white", bold=True, size=9.4, ha="center")
        txt(x + 7.6, 6.4, title, c=DARK, bold=True, size=13)
    # inner dividers between the two levels, inset and dashed
    for xa, xb in [("b", 0), ("c", 0)]:
        x0, _, x1, _ = CARDS[xa]
        line([x0 + 2.4, x1 - 2.4], [69.5, 69.5], "#C3D3D8", .8, style=(0, (3.2, 2.4)))
    # very light flow hint from the molecular anchor into the two scales
    arrow(61.6, 40.0, 65.0, 31.5)
    arrow(61.6, 60.0, 65.0, 83.0)
    # footer band with the central question
    ax.add_patch(Rectangle((3.0, 129.4), 278.0, 10.6, fc=BAND_FC, ec="none", zorder=.4))
    ax.plot([3.0, 281.0], [129.4, 129.4], color=BAND_LINE, lw=.9, zorder=2)
    line([5, 279], [129.4, 129.4], "#B9CDD3", .6, z=2)
    dot(7.0, 135.0, "D", ORANGE, 4.2)
    txt(10.0, 135.0, "Schematic substitution", c="#5E7078", size=9.3)
    txt(281, 135.0, "Do affinity gains carry over to individual mutations?",
        c=DARK, bold=True, ha="right", size=12.2)
elif MODE.startswith("nature"):
    for x, tag, title in PANELS:
        txt(x - 0.6, 6.6, tag, c=DARK, bold=True, size=12.6)
        txt(x + 4.4, 6.6, title, c=DARK, bold=True, size=13.0)
    if MODE == "nature2":
        # hairline grouping rules: header, the two levels and the closing statement
        line([5, 279], [12.6, 12.6], "#DCE3E7", .5)
        for x0, x1 in [(69, 178), (191, 279)]:
            line([x0, x1], [70.6, 70.6], "#E1E7EA", .45)
        line([5, 279], [129.8, 129.8], "#C9D4D9", .6)
    # thin flow arrows: context -> two comparison scales; scales -> strategy detail
    arrow(60.2, 40.0, 66.0, 29.5, "#7E8F96", .8)
    arrow(60.2, 62.0, 66.0, 84.5, "#7E8F96", .8)
    arrow(179.2, 32.0, 187.0, 43.0, "#7E8F96", .8)
    arrow(179.2, 88.0, 187.0, 96.0, "#7E8F96", .8)
    txt(5, 136.0, "Do affinity gains carry over to individual mutations?", c=DARK, bold=True, size=12.6)
    dot(244.5, 136.0, "D", ORANGE, 3.6)
    txt(279.0, 136.0, "Schematic substitution", c="#65757C", size=9.0, ha="right")
elif MODE in ("depth1", "depth3"):
    # 空间纵深：淡背景分区 -> 白色内容卡片（带柔和阴影）-> 前景元素
    for x0, y0, x1, y1 in [(2.6, 16.0, 61.6, 127.4), (64.4, 16.0, 181.6, 127.4), (184.4, 16.0, 281.4, 127.4)]:
        ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                    boxstyle="round,pad=0,rounding_size=2.6",
                                    fc="#F1F6F9", ec="none", zorder=.3))
    def soft_card(x0, y0, x1, y1):
        ax.add_patch(FancyBboxPatch((x0 + .9, y0 + 1.2), x1 - x0, y1 - y0,
                                    boxstyle="round,pad=0,rounding_size=1.5",
                                    fc="#0B1B22", ec="none", alpha=.07, zorder=.6))
        ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                    boxstyle="round,pad=0,rounding_size=1.5",
                                    fc="white", ec="#DDE7EB", lw=.6, zorder=.8))
    if MODE == "depth1":
        soft_card(3.6, 17.2, 60.0, 126.4)
    soft_card(66.2, 17.2, 180.0, 68.9)
    soft_card(66.2, 70.9, 180.0, 126.4)
    soft_card(186.4, 17.2, 279.8, 68.9)
    soft_card(186.4, 70.9, 279.8, 126.4)
    for x, tag, title in PANELS:
        txt(x - 0.8, 6.2, tag, c=DARK, bold=True, size=12.4)
        txt(x + 4.2, 6.2, title, c=DARK, bold=True, size=13.0)
    arrow(60.2, 40.0, 66.0, 29.5, "#7E8F96", .8)
    arrow(60.2, 62.0, 66.0, 84.5, "#7E8F96", .8)
    arrow(179.2, 32.0, 187.0, 43.0, "#7E8F96", .8)
    arrow(179.2, 88.0, 187.0, 96.0, "#7E8F96", .8)
    txt(5, 136.2, "Do affinity gains carry over to individual mutations?", c=DARK, bold=True, size=13.0)
    dot(244.5, 136.2, "D", ORANGE, 3.6)
    txt(279.0, 136.2, "Schematic substitution", c="#65757C", size=9.0, ha="right")
elif MODE == "depth2":
    # 内容层级：两个比较层级用不同深浅的底带体现“层层递进”，其余留白
    ax.add_patch(FancyBboxPatch((2.6, 16.0), 59.0, 111.4, boxstyle="round,pad=0,rounding_size=2.4",
                                fc="#F5F9FB", ec="none", zorder=.3))
    for x0, x1 in [(64.4, 181.6), (184.4, 281.4)]:
        ax.add_patch(FancyBboxPatch((x0, 16.0), x1 - x0, 53.6, boxstyle="round,pad=0,rounding_size=2.4",
                                    fc="#F5F9FB", ec="none", zorder=.3))
        ax.add_patch(FancyBboxPatch((x0, 71.0), x1 - x0, 56.4, boxstyle="round,pad=0,rounding_size=2.4",
                                    fc="#E9F1F6", ec="none", zorder=.3))
    for x, tag, title in PANELS:
        txt(x - 0.8, 6.2, tag, c=DARK, bold=True, size=12.4)
        txt(x + 4.2, 6.2, title, c=DARK, bold=True, size=13.0)
    for x in [64.4 + 3.0, 184.4 + 3.0]:
        txt(x, 88.0, "level 2", c="#7E929B", size=7.2, ha="left")
    for x in [64.4 + 3.0, 184.4 + 3.0]:
        txt(x, 26.5, "level 1", c="#7E929B", size=7.2, ha="left")
    arrow(60.2, 40.0, 66.0, 29.5, "#7E8F96", .8)
    arrow(60.2, 62.0, 66.0, 84.5, "#7E8F96", .8)
    arrow(179.2, 32.0, 187.0, 43.0, "#7E8F96", .8)
    arrow(179.2, 88.0, 187.0, 96.0, "#7E8F96", .8)
    txt(5, 136.2, "Do affinity gains carry over to individual mutations?", c=DARK, bold=True, size=13.0)
    dot(244.5, 136.2, "D", ORANGE, 3.6)
    txt(279.0, 136.2, "Schematic substitution", c="#65757C", size=9.0, ha="right")
elif MODE == "mainsec":
    # 主次层次：白底，放大主标题与关键结论，弱化次要注释，靠字号/灰度拉主次
    for x, tag, title in PANELS:
        txt(x - 0.8, 6.2, tag, c=DARK, bold=True, size=12.2)
        txt(x + 4.2, 6.2, title, c=DARK, bold=True, size=14.4)
    arrow(60.2, 40.0, 66.0, 29.5, "#A6B4BA", .75)
    arrow(60.2, 62.0, 66.0, 84.5, "#A6B4BA", .75)
    arrow(179.2, 32.0, 187.0, 43.0, "#A6B4BA", .75)
    arrow(179.2, 88.0, 187.0, 96.0, "#A6B4BA", .75)
    txt(5, 136.0, "Do affinity gains carry over to individual mutations?", c=DARK, bold=True, size=14.2)
    dot(243.0, 136.0, "D", ORANGE, 4.4)
    txt(279.0, 136.0, "Schematic substitution", c="#8A99A0", size=8.6, ha="right")
elif MODE == "progression":
    # 内容递进：两层比较尺度用“由浅到深”的色阶底带；三个评价层用阶梯色块
    for x0, x1 in [(64.4, 181.6), (184.4, 281.4)]:
        ax.add_patch(FancyBboxPatch((x0, 16.0), x1 - x0, 53.6, boxstyle="round,pad=0,rounding_size=2.4",
                                    fc="#F4F8FB", ec="none", zorder=.3))
        ax.add_patch(FancyBboxPatch((x0, 71.0), x1 - x0, 56.4, boxstyle="round,pad=0,rounding_size=2.4",
                                    fc="#E2EDF4", ec="none", zorder=.3))
    ax.add_patch(FancyBboxPatch((2.6, 16.0), 59.0, 111.4, boxstyle="round,pad=0,rounding_size=2.4",
                                fc="#F7FAFC", ec="none", zorder=.3))
    for x, tag, title in PANELS:
        txt(x - 0.8, 6.2, tag, c=DARK, bold=True, size=12.4)
        txt(x + 4.2, 6.2, title, c=DARK, bold=True, size=13.0)
    arrow(60.2, 40.0, 66.0, 29.5, "#7E8F96", .8)
    arrow(60.2, 62.0, 66.0, 84.5, "#7E8F96", .8)
    arrow(179.2, 32.0, 187.0, 43.0, "#7E8F96", .8)
    arrow(179.2, 88.0, 187.0, 96.0, "#7E8F96", .8)
    txt(5, 136.2, "Do affinity gains carry over to individual mutations?", c=DARK, bold=True, size=13.0)
    dot(244.5, 136.2, "D", ORANGE, 3.6)
    txt(279.0, 136.2, "Schematic substitution", c="#65757C", size=9.0, ha="right")
elif MODE == "glossy":
    # 立体/设计感：渐变底板 + 悬浮白卡（双层柔和投影 + 顶部高光）
    if not FLAT:
        from matplotlib.colors import LinearSegmentedColormap
        zone_cmap = LinearSegmentedColormap.from_list("zone", ["#FBFDFE", "#E8F1F6"])
        grad = np.linspace(0, 1, 64).reshape(-1, 1)
        for x0, y0, x1, y1 in [(2.6, 16.0, 61.6, 127.4), (64.4, 16.0, 181.6, 127.4), (184.4, 16.0, 281.4, 127.4)]:
            ax.imshow(grad, extent=[x0, x1, y0, y1], cmap=zone_cmap, aspect="auto",
                      interpolation="bilinear", zorder=.2)
        def float_card(x0, y0, x1, y1):
            for dx, dy, al in [(2.4, 2.8, .04), (1.3, 1.6, .06)]:
                ax.add_patch(FancyBboxPatch((x0 + dx, y0 + dy), x1 - x0, y1 - y0,
                                            boxstyle="round,pad=0,rounding_size=2.0",
                                            fc="#0B1B22", ec="none", alpha=al, zorder=.6))
            ax.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                        boxstyle="round,pad=0,rounding_size=2.0",
                                        fc="white", ec="#E3EBEF", lw=.6, zorder=.8))
            ax.plot([x0 + 3.0, x1 - 3.0], [y0 + .9, y0 + .9], color="white", lw=.9, zorder=.9, alpha=.9)
        float_card(2.0, 16.8, 60.6, 126.8)
        float_card(66.0, 16.8, 180.8, 68.9)
        float_card(66.0, 71.0, 180.8, 126.8)
        float_card(188.0, 16.8, 280.2, 68.9)
        float_card(188.0, 71.0, 280.2, 126.8)
    for x, tag, title in PANELS:
        txt(x - 0.6, 6.4, tag, c=DARK, bold=True, size=12.4)
        txt(x + 4.4, 6.4, title, c=DARK, bold=True, size=13.2)
    arrow(60.4, 40.0, 65.8, 29.5, "#7E8F96", .8)
    arrow(60.4, 62.0, 65.8, 84.5, "#7E8F96", .8)
    arrow(179.2, 32.0, 186.8, 43.0, "#7E8F96", .8)
    arrow(179.2, 88.0, 186.8, 96.0, "#7E8F96", .8)
    txt(5, 136.2, "Do affinity gains carry over to individual mutations?", c=DARK, bold=True, size=13.0)
    dot(244.5, 136.2, "D", ORANGE, 3.6)
    txt(279.0, 136.2, "Schematic substitution", c="#65757C", size=9.0, ha="right")
else:
    for x, tag, title in PANELS:
        ax.add_patch(Rectangle((x - 1.4, 2.8), 4.9, 4.9, fc="white", ec=BLUE, lw=.85, zorder=6))
        txt(x + 1.05, 5.35, tag, c=BLUE, bold=True, size=8.9, ha="center")
        txt(x + 6.2, 6.0, title, c=DARK, bold=True, size=13)
    line([5, 279], [13.0, 13.0], "#7E979D", 1.0)
    for x in [63.5, 184.5]:
        line([x, x], [15.4, 126.8], "#A7BBC1", 1.0)
    for x0, x1 in [(69, 178), (191, 279)]:
        line([x0, x1], [69.5, 69.5], "#B9CBD0", .85)
    line([5, 279], [129.0, 129.0], "#16323C", 1.4)
    dot(7.0, 135.5, "D", ORANGE, 4.2)
    txt(10.0, 135.5, "Schematic substitution", c="#5E7078", size=9.3)
    txt(279, 135.5, "Do affinity gains carry over to individual mutations?",
        c=DARK, bold=True, ha="right", size=12.2)

# ---------------------------------------------------------------- figure content
def pair_2c4r(m):
    """The 2C4R Phe57 CZ - G10 O2' contact discussed in the Results."""
    f57 = [r for r in m["L"] if r.id[1] == 57][0]
    g10 = [r for r in m["R"] if r.id[1] == 10][0]
    return (2.80, f57, g10, f57["CZ"], g10["O2'"])


struct_img("2C4R_ribbon.png", 31, 38, 36, clip="a", bleed=0.85)
txt(31, 62, "2C4R · RNase E–ssRNA", c=DARK, ha="center", bold=True, size=10.2)
struct_img("2C4R_contact.png", 31, 92, 36, clip="a", bleed=0.85)
txt(31, 118.5, "Phe57–G10 contact · 2.80 Å", c=GRAY, ha="center", size=9.4)

cap(69, 21, "Across binding systems")
txt(176.4, 21, "Absolute affinity", c=GRAY, ha="right", size=9.6)
for pdb, x, rna, i in [("4QVC", 84, "Hfq–OxyS", 1), ("2ZZN", 123, "Trm5–tRNA", 2), ("4NGD", 162, "Dicer–siRNA", 3)]:
    struct_img(f"{pdb}_ribbon.png", x, 33.5, 21, clip="b1", bleed=0.62)
    txt(x, 47.2, f"{pdb} · {rna}", c=GRAY, ha="center", size=10.0)
    line([x, x], [49.8, 52.3], "#A0B5B8", .65)
    name(x, 55.8, rf"$\mathrm{{p}}K_{{d,{i}}}$", DARK)
line([84, 162], [61, 61], "#A0B5B8", .8)
for x in [84, 123, 162]:
    line([x, x], [58.2, 61], "#A0B5B8", .65)
txt(176.4, 65.4, "Pooled endpoint error", c=DARK, bold=True, ha="right", size=10)

cap(69, 75.2, "Within one matched system")
txt(176.4, 80.6, "Same RNA and assay", c=GRAY, ha="right", size=9.6)
labels = ["WT", "F57A", "R373A", "D346N", "N305L"]
positions = [78, 101, 124, 147, 170]
for i, x in enumerate(positions):
    schematic(x, 90, .78, 0 if i == 0 else (i % 4) + 1)
    txt(x, 98, labels[i], c=DARK, ha="center", size=9.4)
    txt(x, 103.4, rf"$\mathrm{{p}}K_{{d,\mathrm{{{labels[i]}}}}}$", c=DARK, ha="center", size=9.0)
    line([x, x], [105.8, 109], "#88A3AA", .7)
line([78, 170], [109, 109], "#88A3AA", .85)
for i, x in enumerate(positions[1:], 1):
    ax.add_patch(Circle((x, 109), 1.6, fc="white", ec=TEAL, lw=.8, zorder=6))
    name(x, 108.94, "−", DARK)
    line([x, x], [110.6, 112.2], "#88A3AA", .65)
    txt(x, 115.5, rf"$\Delta\mathrm{{p}}K_{{d,\mathrm{{{labels[i]}}}}}$", c=DARK, ha="center", size=8.6)
for x, label in [(104, "Mean"), (134, "Order"), (164, "Spread")]:
    txt(x, 124, label, c=DARK, ha="center", size=10.5)
line([99, 109], [120.7, 120.7], TEAL, .85)
for x in [100, 104, 108]:
    dot(x, 120.7, "o", TEAL, 2.3)
for x, y in [(130, 121.6), (134, 120.7), (138, 119.8)]:
    dot(x, y, "o", TEAL, 2.3)
line([159, 169], [120.7, 120.7], TEAL, .85)
for x in [159, 169]:
    line([x, x], [119.7, 121.7], TEAL, .85)

cap(191, 21, "Compare matched predictions")
token_row(198, 31, "AKGR")
token_row(198, 35.2, "ACGU")
graph(233, 33)
art_c = ax.imshow(ribbon_image(), extent=[260, 272, 40, 26], aspect="auto", zorder=4)
token_row(260, 42.2, "AKGR")
for x, title in [(204, "Sequence"), (233, "Structure"), (267, "Multimodal")]:
    name(x, 49, title)
    line([x, x], [51.5, 54.5], "#A0B5B8", .65)
line([204, 267], [54.5, 54.5], "#A0B5B8", .8)
txt(233.6, 59.2, "Same-member populations", c=DARK, ha="center", size=10.5)
for x, title in [(203, "Absolute"), (234, "Mutation"), (266, "Centered")]:
    txt(x, 65.4, title, c=GRAY, ha="center", size=9.6)

cap(191, 75.2, "Locate the mutation response")
txt(191, 84.5, "PRA-Pred", c=GRAY, size=9.6)
txt(216, 84.5, "WT", c=TEAL, size=9.4)
txt(225.5, 84.5, "Mut", c=ORANGE, size=9.4)
txt(232, 84.5, r"$\widehat{\Delta\Delta G}=\sum_j\beta_j\Delta x_j$", c=DARK, size=9.8)

for x, label in [(206, "Atom"), (232, "Residue"), (258, "Image")]:
    if label == "Atom":
        graph(x, 98)
    elif label == "Residue":
        curve([(x - 5, 98), (x - 3, 94), (x, 102), (x + 5, 97)], TEAL, .9)
        for dx in [-3, 0, 3]:
            dot(x + dx, 98, "o", TEAL, 2.3, fill=False)
    else:
        for dx, dy in [(-2.8, -1.2), (-1.4, 0), (0, 1.2)]:
            ax.add_patch(Rectangle((x + dx, 96 + dy), 6, 3.5, fc="#EDF3F4",
                                  ec="#7F9FA7", lw=.5, zorder=4))
    txt(x, 104.5, label, c=GRAY, ha="center", size=9.4)
    for yy, c in [(108, TEAL), (109.7, ORANGE)]:
        for dx in [-3, -1, 1, 3]:
            line([x + dx, x + dx], [yy - .35, yy + .35], c, 1.0, z=4)
    line([x, x], [111.3, 114.3], "#95AFB5", .65)
line([206, 269], [114.3, 114.3], "#95AFB5", .8)
txt(274.0, 114.3, "Fusion", c=DARK, ha="center", size=9.4)
txt(233.2, 123.5, "M2-PRNet · input / encoder / readout", c=DARK, ha="center", size=10)

export("molecular_context_grounded")
print("rendered", "flat" if FLAT else MODE, "->", DEST)
