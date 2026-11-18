"""Render Supplementary Figure E1 from recomputed seed and family statistics."""
import numpy as np

NAMES = ("sf3a1", "sepsecs", "npl3")


def make_figure(primary, seed_rows, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D
    plt.rcParams.update({'font.family':'Arial', 'font.size':8, 'axes.titlesize':10,
                         'axes.labelsize':8, 'pdf.fonttype':42, 'ps.fonttype':42,
                         'svg.fonttype':'none', 'axes.spines.top':False, 'axes.spines.right':False})
    fig, axes = plt.subplots(1, 3, figsize=(7.0, 3.6), sharey=True)
    fig.subplots_adjust(left=.17, right=.985, top=.74, bottom=.32, wspace=.24)
    keys = ['absolute_MAE', 'ordinary_MAE', 'cMAE']
    titles = ['Absolute affinity', 'Mutation effect', 'Local variation']
    y = [3, 2, 1, 0]
    for a, key, title, letter in zip(axes, keys, titles, 'abc'):
        a.axvline(0, color='#999999', linestyle='--', linewidth=.8, zorder=0)
        for i, name in enumerate(NAMES):
            by_method = {m:{r['seed']:r[key] for r in seed_rows if r['family']==name and r['method']==m} for m in ('copra','baseline')}
            values = [by_method['baseline'][s]-by_method['copra'][s] for s in sorted(by_method['copra'])]
            assert all(-.12 <= v <= .42 for v in values), (name, key, values)
            a.scatter(values, np.array([y[i]]*3)+[-.09,0,.09], s=16, facecolors='white', edgecolors='#89969D', linewidths=.8, zorder=2)
            mean = primary['contrasts'][key]['family_values'][name]
            a.scatter([mean], [y[i]], s=28, c='#237A8E', zorder=3)
        contrast = primary['contrasts'][key]
        lo, hi = contrast['bootstrap95']
        mid = contrast['mean']
        a.errorbar(mid, 0, xerr=[[mid-lo],[hi-mid]], fmt='D', markersize=5,
                   color='#C86A43', linewidth=1.2, capsize=3, zorder=4)
        a.axhline(.5, color='#E1E4E6', linewidth=.7)
        a.set_title(title, loc='left', pad=10, fontweight='bold')
        a.text(-.08, 1.14, letter, transform=a.transAxes, fontweight='bold', fontsize=12)
        a.set_xlim(-.12,.42)
        a.set_xticks([-.1,0,.1,.2,.3,.4])
        a.set_ylim(-.55,3.55)
        a.tick_params(labelsize=7.5)
        a.set_xlabel('Error reduction (pKd)')
    axes[0].set_yticks(y, ['SF3A1 (7 pairs)', 'SEPSECS (11)', 'Npl3 (6)', 'Family mean'])
    fig.text(.17,.96,'Three source-matched external families',fontsize=12,fontweight='bold',color='#202A30')
    fig.text(.17,.88,'24 quantitative pairs; fixed weights; additive minus controlled CoPRA',fontsize=8.5,color='#607078')
    handles=[Line2D([],[],marker='o',linestyle='',color='#89969D',markerfacecolor='white',markersize=4,label='Seed'),
             Line2D([],[],marker='o',linestyle='',color='#237A8E',markersize=5,label='Family mean across seeds'),
             Line2D([],[],marker='D',linestyle='-',color='#C86A43',markersize=4,label='Three-family mean and 95% CI')]
    fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.53,.09),ncol=3,frameon=False,fontsize=7.5,handlelength=1.3,columnspacing=1.5)
    fig.text(.17,.035,'Positive values favor CoPRA. CI: 10,000 paired family-bootstrap draws (3 families).',fontsize=7.5,color='#607078')
    fig.canvas.draw()
    renderer=fig.canvas.get_renderer()
    clipped=[]
    for text in fig.findobj(matplotlib.text.Text):
        if not text.get_visible() or not text.get_text():continue
        box=text.get_window_extent(renderer)
        if box.x0 < -1 or box.y0 < -1 or box.x1 > fig.bbox.x1+1 or box.y1 > fig.bbox.y1+1:
            clipped.append(text.get_text())
    assert not clipped, clipped
    fig.savefig(output, dpi=300)
    plt.close(fig)
