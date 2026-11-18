"""Figure 2, wide manuscript reflow: 178.0 mm canvas, native typography preserved.

Native design: 142.4 x 223 mm, rendered by
outputs/figure2_muted_house_palette_20261102_v2_restyled/render_figure2.py.
This renderer keeps every text string, every native point size (7.5/8/9/11 pt),
the frozen plot_data.json values, all axis limits/ticks, the palette and the
panel letters A-D.  Only the layout coordinates are recomposed:

  row 1 : A (error decomposition, left)  +  B (4x4 residual corner plot, right)
  row 2 : C (52-family x seed bar + atlas) spanning the full text column
  row 3 : D (joint bootstrap + marginal densities + forest) full width

Each panel keeps its native internal coordinate system; a small affine map per
panel converts native mm into the new 178 mm-wide canvas.  Because the maps only
translate/scale positions (never text), the native matplotlib point sizes are
untouched.
"""
from pathlib import Path
import json
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm, LinearSegmentedColormap
from matplotlib.ticker import FuncFormatter
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import gaussian_kde

O = Path(__file__).resolve().parent
D = json.loads((O / 'plot_data.json').read_text(encoding='utf-8'))
v5_panels = '--v5-panels' in sys.argv
W, H = 178.0, 199.0
from bright_palette import BLUE, TEAL, CORAL, AMBER, INK, SECONDARY, AXIS, GRID, LIGHT_BLUE, LIGHT_TEAL
GRAY = SECONDARY
plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'axes.labelsize': 8,
    'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5, 'axes.linewidth': .55,
    'xtick.major.width': .5, 'ytick.major.width': .5, 'xtick.major.size': 2,
    'ytick.major.size': 2, 'xtick.minor.visible': False, 'ytick.minor.visible': False,
    'text.color': INK, 'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK,
    'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none',
    'mathtext.fontset': 'dejavusans'})
fig = plt.figure(figsize=(W / 25.4, H / 25.4), facecolor='white')


# --------------------------------------------------------------------------
# absolute-mm drawing helpers (the new canvas is native mm too)
# --------------------------------------------------------------------------
def AX(x, y, w, h):
    """Axes from absolute mm (x right, y down from the page top)."""
    return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])


def TX(x, y, s, **kw):
    return fig.text(x / W, 1 - y / H, s, va=kw.pop('va', 'top'), **kw)


def TITLE(x, y, letter, label):
    TX(x, y, letter, fontsize=11, fontweight='bold')
    TX(x + 5.2, y + .65, label, fontsize=9, fontweight='bold')


def LEG(handles, x, y, **kw):
    """Upper-left anchored figure legend at absolute mm."""
    return fig.legend(handles=handles, loc='upper left', bbox_to_anchor=(x / W, 1 - y / H),
                      frameon=False, fontsize=7.5, borderaxespad=0, **kw)


def style(ax, box=False):
    ax.tick_params(pad=1.5)
    for s in ax.spines.values():
        s.set_color(AXIS); s.set_linewidth(.5)
    if not box:
        ax.spines[['top', 'right']].set_visible(False)


def hdr_levels(z):
    v = np.sort(z.ravel())[::-1]
    cs = np.cumsum(v) / v.sum()
    levels = [v[min(np.searchsorted(cs, p), len(v) - 1)] for p in [.95, .8, .5]]
    return np.unique([*levels, z.max() * 1.0001])


# --------------------------------------------------------------------------
# per-panel native-mm -> canvas-mm maps
#   row 1 : A native content column x=11, title x=3, y origin 4
#           B native label/title anchor x=60, grid x0=68, y origin 4
#   row 2 : C native bar top-left (19, 111); title anchor y=100
#   row 3 : D native content origin x=19, title anchor y=166
# --------------------------------------------------------------------------
R1_Y = 2.0
# A keeps its native absolute-mm stack.  The native geo axes box held a
# 25.33 mm-tall equal-aspect data box inside a 34 mm slot; in the wide canvas
# the plot is taller, so the fixed-height text rows below it (x-label, MSE
# formula, legends, MSE bar) are positioned at their NATIVE millimetre
# clearances from the bottom of the equal-aspect data box instead of being
# squashed with the axes.  That keeps the formula clear of the x-axis label.
A_X, A_SX = 12.0, 1.42                      # A content column, 38 mm native -> 54 mm
A_SY = 0.91                                 # geo / bar height scale
A_GEO_Y = R1_Y + 19 * A_SY                  # native geo slot top y=23 (19 below title)
A_GEO_H = 34 * A_SY
A_GEO_BOT = A_GEO_Y + A_GEO_H
A_LEG1_Y = A_GEO_Y - 6.0                    # native legend anchor y=17 is 6 above slot top
A_EQ_Y = A_GEO_BOT + 9.04                   # formula: native geo-box bottom + 9.04
A_LEG2_Y = A_GEO_BOT + 17.84                # shared/difference legend (native y=70.5)
A_PART_Y = A_GEO_BOT + 23.34                # MSE bar top (native y=76)
A_PART_H = 13 * A_SY
B_XB, B_SBX, B_SY = 79.36, 1.2944, 0.995    # B anchor x=60; grid ends at 172 mm
C_TY, C_X, C_SX, C_SY = 95.0, 22.0, 1.31579, 0.84
C_Y = C_TY + 11 * C_SY                      # C native bar top y=111
D_TY, D_X, D_SX, D_SY = 148.3, 19.0, 1.34211, 0.90
D_HEAD2_X = 110.7                           # "Mean reduction" header sits over the forest cluster


def mapB(x, y, w=0.0, h=0.0):
    return B_XB + (x - 60) * B_SBX, R1_Y + (y - 4) * B_SY, w * B_SBX, h * B_SY


def mapC(x, y, w=0.0, h=0.0):
    return C_X + (x - 19) * C_SX, C_Y + (y - 111) * C_SY, w * C_SX, h * C_SY


def mapD(x, y, w=0.0, h=0.0):
    return D_X + (x - 19) * D_SX, D_TY + (y - 166) * D_SY, w * D_SX, h * D_SY


def axB(x, y, w, h):
    return AX(*mapB(x, y, w, h))


def axC(x, y, w, h):
    return AX(*mapC(x, y, w, h))


def axD(x, y, w, h):
    return AX(*mapD(x, y, w, h))


def txB(x, y, s, **kw):
    xx, yy, _, _ = mapB(x, y)
    return TX(xx, yy, s, **kw)


def txC(x, y, s, **kw):
    xx, yy, _, _ = mapC(x, y)
    return TX(xx, yy, s, **kw)


def txD(x, y, s, **kw):
    xx, yy, _, _ = mapD(x, y)
    return TX(xx, yy, s, **kw)


# ==========================================================================
# A: actual residuals in decomposed coordinates plus mean composition.
# ==========================================================================
TITLE(3, R1_Y, 'A', 'Error decomposition')
geo = AX(A_X, A_GEO_Y, 38 * A_SX, A_GEO_H)
rot = D['rotated_residuals']; rx, ry = np.asarray(rot['x']), np.asarray(rot['y'])
for model in ['copra', 'baseline']:
    density = np.asarray(rot['densities'][model]); lev = hdr_levels(density)
    if model == 'copra':
        geo.contourf(rx, ry, density, levels=lev, colors=['#E8EFF6', '#C3D6E7', '#8FB2D0'])
    geo.contour(rx, ry, density, levels=lev[:-1], colors=BLUE if model == 'copra' else TEAL,
                linestyles='solid' if model == 'copra' else 'dashed', linewidths=.55)
for radius in [1, 2, 3]:
    theta = np.linspace(0, 2 * np.pi, 240)
    geo.plot(radius * np.cos(theta), radius * np.sin(theta), color=GRID, lw=.4, ls=':', zorder=0)
    geo.text(radius * .5 + .15, -radius * .866, str(radius ** 2), fontsize=7.5, color=SECONDARY, va='center',
             bbox=dict(facecolor='white', edgecolor='none', alpha=.9, pad=.15))
geo.axhline(0, color=GRID, lw=.35, zorder=0); geo.axvline(0, color=GRID, lw=.35, zorder=0)
geo.set(xlim=(-6, 6), ylim=(-4, 4), xticks=[-6, 0, 6], yticks=[-3, 0, 3], aspect='equal')
geo.set_xlabel(r'$c$ (pKd)', labelpad=0)
geo.set_ylabel(r'$q/2$ (pKd)', labelpad=1)
style(geo)
LEG([Line2D([], [], color=TEAL, ls='--', lw=.7, label='Baseline'),
     Line2D([], [], color=BLUE, lw=.8, label='CoPRA')],
    A_X, A_LEG1_Y, ncol=2, handlelength=1, handletextpad=.4, columnspacing=.8)
TX(A_X, A_EQ_Y, r'$\mathrm{MSE}_{\mathrm{endpoint}}=\langle c^2\rangle+\langle q^2/4\rangle$', fontsize=8)
LEG([Patch(facecolor=BLUE, label='Shared'), Patch(facecolor=TEAL, edgecolor=BLUE, hatch='///', label='Difference')],
    A_X, A_LEG2_Y, ncol=2, handlelength=.9, handletextpad=.4, columnspacing=.8)
part = AX(A_X, A_PART_Y, 38 * A_SX, A_PART_H)
for y, (m, label) in zip([1, 0], [('baseline', 'Baseline'), ('copra', 'CoPRA')]):
    shared, diff, total = D['means'][m]
    part.barh(y, shared, height=.58, color=BLUE, linewidth=0)
    part.barh(y, diff, left=shared, height=.58, color=TEAL, hatch='///', edgecolor=BLUE, linewidth=.3)
    part.text(.08, y, label, va='center', fontsize=7.5, color=INK)
    part.text(total + .08, y, f'{total:.3f}', va='center', fontsize=7.5)
part.set(xlim=(0, 3.5), ylim=(-.6, 1.6), yticks=[], xticks=[0, 1, 2, 3])
part.tick_params(axis='y', length=0, pad=3)
part.set_xlabel(r'MSE (pKd$^2$)', fontsize=7.5, labelpad=1)
style(part)

# ==========================================================================
# B: four actual residual variables, not replicate measurements.
# ==========================================================================
blue_colors = ['#E8EFF6', '#C3D6E7', '#8FB2D0']
if not v5_panels:
    TITLE(B_XB, R1_Y, 'B', 'Residual dependence')
    txB(103, 6, '338 pairs · 52 families', fontsize=7.5)
    variables = ['Baseline\nWT', 'Baseline\nMut', 'CoPRA\nWT', 'CoPRA\nMut']
    grid = np.asarray(D['pairplot']['grid'])
    hist_edges = np.asarray(D['pairplot']['hist_edges'])
    histmax = np.max(D['pairplot']['histograms']) * 1.12
    cell, gap, x0, y0 = 15.0, 1.2, 68, 21
    teal_colors = ['#F8EFE2', '#EAD3B3', '#D5A86F']
    for i in range(4):
        txB(x0 - 8, y0 + i * (cell + gap) + cell / 2, variables[i], ha='right', va='center', fontsize=7.5)
        txB(x0 + i * (cell + gap) + cell / 2, 12, variables[i], ha='center', fontsize=7.5)
        for j in range(4):
            ax = axB(x0 + j * (cell + gap), y0 + i * (cell + gap), cell, cell)
            if i == j:
                hist = np.asarray(D['pairplot']['histograms'][i])
                ax.stairs(hist, hist_edges, fill=True, color=BLUE if i >= 2 else TEAL, linewidth=.5)
                ax.set(ylim=(0, histmax), yticks=[])
            else:
                entry = D['pairplot']['pairs'][f'{max(i, j)},{min(i, j)}']
                density = np.asarray(entry['density'])
                if i < j: density = density.T
                levels = hdr_levels(density)
                ax.contourf(grid, grid, density, levels=levels, colors=blue_colors if i > j else teal_colors, antialiased=True)
                ax.contour(grid, grid, density, levels=levels[:-1], colors=BLUE if i > j else TEAL, linewidths=.3, alpha=.65)
                ax.plot([-7, 7], [-7, 7], color=AXIS, lw=.45, ls=(0, (2, 2)), zorder=0)
                ax.set(ylim=(-7, 7), yticks=[-6, 0, 6])
                if i > j:
                    ax.text(.055, .94, r'$r_w=' + f"{entry['r_weighted']:.2f}" + '$', transform=ax.transAxes, va='top', fontsize=7.5,
                            bbox=dict(facecolor='white', edgecolor='none', alpha=.88, pad=.3))
            ax.set(xlim=(-7, 7), xticks=[-6, 0, 6])
            if i != 3: ax.set_xticklabels([])
            if j != 0 or i == j: ax.set_yticklabels([])
            style(ax)

    txB(100, 92, 'Endpoint residuals (pKd)', ha='center', fontsize=7.5)

# ==========================================================================
    # C: total family gain followed by the full 52-family x 3-seed component atlas.
# ==========================================================================
    TITLE(3, C_TY, 'C', 'Family × seed decomposition')
    txC(85, 101, 'Improved: 31/52 shared · 41/52 diff.', fontsize=7.5)
    order = np.asarray(D['atlas']['order'])
    gain = np.asarray(D['atlas']['seed_reductions'])[order]
    mean_gain = gain.mean(axis=1)
    total = mean_gain[:, 2]
    bar = axC(19, 111, 114, 13)
    bar.bar(np.arange(52), total, width=.77, color=np.where(total >= 0, BLUE, CORAL), linewidth=0)
    bar.axhline(0, color=GRAY, lw=.5)
    bar.set(xlim=(-.5, 51.5), ylim=(-12, 12), yticks=[-10, 0, 10], xticks=[])
    bar.set_ylabel('Total\ngain', fontsize=7.5, labelpad=2)
    style(bar)
    atlas = axC(19, 127, 114, 21)
    matrix = np.vstack([gain[:, :, 0].T, gain[:, :, 1].T])
    vmax = float(np.ceil(np.abs(matrix).max()))
    norm = SymLogNorm(linthresh=.03, linscale=.6, vmin=-vmax, vmax=vmax, base=10)
    cmap = LinearSegmentedColormap.from_list('gain', [CORAL, '#D9A08F', '#FFFFFF', '#B3CCDF', BLUE])
    im = atlas.pcolormesh(np.arange(53) - .5, np.arange(7) - .5, matrix, cmap=cmap, norm=norm, edgecolors='white', linewidth=.15)
    atlas.set(xlim=(-.5, 51.5), ylim=(5.5, -.5), yticks=np.arange(6), yticklabels=['S1', 'S2', 'S3'] * 2,
              xticks=[0, 12, 25, 38, 51], xticklabels=['1', '13', '26', '39', '52'])
    atlas.tick_params(axis='y', length=0, pad=3)
    atlas.tick_params(axis='x', length=2, pad=1.5)
    atlas.axhline(2.5, color=INK, lw=.65)
    for s in atlas.spines.values(): s.set_visible(False)
    TX(2.0, C_Y + (132.25 - 111) * C_SY, 'Shared', fontsize=7.5, ha='left', va='center')
    TX(2.0, C_Y + (142.75 - 111) * C_SY, 'Difference', fontsize=7.5, ha='left', va='center')
    txC(19, 155, 'Family rank by total MSE reduction', fontsize=7.5)
    ca = axC(104, 155, 29, 1.8)
    cb = fig.colorbar(im, cax=ca, orientation='horizontal', ticks=[-10, 0, 10])
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=7.5, length=1.5, pad=1)
    cb.ax.xaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'{v:g}'))
    cb.ax.minorticks_off()
    if cb.solids is not None: cb.solids.set_rasterized(False)
    txC(118.5, 151, 'Gain (pKd²)', ha='center', fontsize=7.5)


# ==========================================================================
# D: dependence between bootstrap estimates, alongside marginal intervals.
# ==========================================================================
TITLE(3, D_TY, 'D', 'Joint uncertainty and mean gains')
txD(19, 175, 'Shared–difference coupling', fontsize=7.5)
TX(D_HEAD2_X, mapD(83, 175)[1], 'Mean reduction · 95% CI', fontsize=7.5)
b = D['joint_bootstrap']; bx, by = np.asarray(b['x']), np.asarray(b['y']); bd = np.asarray(b['density'])
joint = axD(19, 187, 46, 24)
lev = hdr_levels(bd)
joint.contourf(bx, by, bd, levels=lev, colors=blue_colors)
joint.contour(bx, by, bd, levels=lev[:-1], colors=BLUE, linewidths=.5)
joint.axvline(0, color=GRAY, lw=.5)
joint.plot(-by, by, color=INK, lw=.6, ls=(0, (3, 2)))
if v5_panels:
    point = [D['effects'][key]['estimate'] for key in ('shared_mse', 'difference_mse')]
else:
    point = mean_gain.mean(axis=0)
joint.scatter(point[0], point[1], marker='D', s=23, color=AMBER, edgecolors=INK, linewidths=.5, zorder=4)
joint.set(xlim=(bx[0], bx[-1]), ylim=(by[0], by[-1]), xticks=[-.5, 0, .5, 1], yticks=[0, .05, .1])
joint.yaxis.set_major_formatter(FuncFormatter(lambda v, pos: f'{v:g}'))
joint.set_xlabel(r'Shared reduction (pKd$^2$)', labelpad=2)
joint.set_ylabel(r'Difference reduction (pKd$^2$)', labelpad=2)
style(joint)
mt = axD(19, 182, 46, 3.5); mt.fill_between(bx, bd.sum(axis=0), color=BLUE, alpha=.5, lw=0)
mt.set(xlim=(bx[0], bx[-1]), ylim=(0, bd.sum(axis=0).max() * 1.1)); mt.axis('off')
mr = axD(66.5, 187, 3.5, 24); mr.fill_betweenx(by, bd.sum(axis=1), color=TEAL, alpha=.5, lw=0)
mr.set(ylim=(by[0], by[-1]), xlim=(0, bd.sum(axis=1).max() * 1.1)); mr.axis('off')
forest = axD(100, 183, 33, 28)
xs = np.asarray(D['effect_density_grid']) if v5_panels else np.linspace(-.8, 1.5, 600)
for y, (key, label, color) in zip([2, 1, 0], [('shared_mse', 'Shared', BLUE), ('difference_mse', 'Difference', TEAL), ('endpoint_mse', 'Total', CORAL)]):
    v = D['effects'][key]
    z = np.asarray(v['density']) if v5_panels else gaussian_kde(v['bootstrap'])(xs)
    z = .35 * z / z.max()
    forest.fill_between(xs, y + .045, y + .045 + z, color=color, alpha=.5, lw=0)
    lo, hi = v['ci95']; estimate = v['estimate']
    forest.plot([lo, hi], [y, y], color=color, lw=1.2)
    forest.plot([lo, lo], [y - .07, y + .07], color=color, lw=.6)
    forest.plot([hi, hi], [y - .07, y + .07], color=color, lw=.6)
    forest.scatter(estimate, y, marker='D', s=10, color=AMBER, edgecolors=INK, linewidths=.4, zorder=4)
    forest.text(estimate, y - .18, f'{estimate:.3f}', ha='center', va='top', fontsize=7.5)
forest.axvline(0, color=GRAY, ls=(0, (3, 2)), lw=.6)
forest.set(xlim=(-.7, 1.35), ylim=(-.5, 2.6), xticks=[-.5, 0, .5, 1], yticks=[2, 1, 0], yticklabels=['Shared', 'Difference', 'Total'])
forest.tick_params(axis='y', length=0, pad=3)
forest.set_xlabel(r'MSE reduction (pKd$^2$)', labelpad=2)
style(forest); forest.spines['left'].set_visible(False)

# ==========================================================================
# visual QA (deterministic layer) + exports
# ==========================================================================
fig.canvas.draw()
sys.path.insert(0,str(Path(__file__).resolve().parent))
from visual_qa import audit_layout
issues = audit_layout(fig)
(O / 'layout_audit.json').write_text(json.dumps({'issues': issues, 'size_mm': [W, H], 'dpi': 400}, indent=2))
assert not issues, issues
fig.savefig(O / 'figure2.png', dpi=400, facecolor='white')
fig.savefig(O / 'figure2.pdf', facecolor='white')
fig.savefig(O / 'figure2.svg', facecolor='white')
fig.savefig(O / 'preview.png', dpi=150, facecolor='white')
print(json.dumps({'status': 'rendered', 'layout_issues': issues, 'size_mm': [W, H]}))
