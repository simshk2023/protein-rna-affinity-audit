"""Render RFdiffusion-style rainbow cartoons with the project's existing PyMOL env."""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PDB = ROOT / 'outputs' / 'length_research_handoff_20260926' / 'DS04_work' / 'sources' / 'PRA310_official' / 'PRA310' / 'PDBs'
OUT = ROOT / 'outputs' / 'problem_framing_ribbon_20261107' / 'renders'
OUT.mkdir(parents=True, exist_ok=True)

RNA_BLUE = [0.00, 0.32, 0.50]


def setup_common():
    from pymol import cmd
    cmd.reinitialize()
    cmd.bg_color('white')
    cmd.set('ray_opaque_background', 0)
    cmd.set('antialias', 2)
    cmd.set('cartoon_fancy_helices', 1)
    cmd.set('cartoon_smooth_loops', 1)
    cmd.set('cartoon_transparency', 0.0)
    cmd.set('ray_shadows', 0)
    cmd.set('ray_trace_mode', 0)
    cmd.set_color('rna_blue', RNA_BLUE)


def render(name, pdb, prot_chains, rna_chains, size=1800, chainwise=False, zoom_extra=2.0,
           rna_sticks=False, transparency=0.0, hide_chains=()):
    from pymol import cmd
    setup_common()
    cmd.load(str(pdb), name)
    cmd.hide('everything')
    if hide_chains:
        cmd.hide('everything', ' or '.join(f'chain {c}' for c in hide_chains))
    cmd.show('cartoon', 'polymer')
    cmd.show('cartoon', 'polymer and name P')
    cmd.color('grey70', 'polymer')
    if chainwise:
        for i, ch in enumerate(prot_chains):
            cmd.spectrum('resi', 'rainbow', f'chain {ch} and polymer')
    else:
        cmd.spectrum('resi', 'rainbow', 'polymer and (' + ' or '.join(f'chain {c}' for c in prot_chains) + ')')
    for ch in rna_chains:
        cmd.show('cartoon', f'chain {ch} and polymer')
        cmd.color('rna_blue', f'chain {ch} and polymer')
        cmd.set('cartoon_transparency', 0.0, f'chain {ch} and polymer')
        if rna_sticks:
            cmd.show('sticks', f'chain {ch} and polymer')
            cmd.set('stick_radius', 0.18, f'chain {ch} and polymer')
    if transparency:
        cmd.set('cartoon_transparency', transparency, 'polymer')
    cmd.orient('polymer')
    cmd.zoom('polymer', zoom_extra)
    cmd.ray(size, size)
    out = OUT / f'{name}_ribbon.png'
    cmd.png(str(out), dpi=300)
    print('wrote', out, out.stat().st_size)
    return out


def render_contact(size_w=1800, size_h=1500):
    from pymol import cmd
    setup_common()
    cmd.load(str(PDB / '2C4R.pdb'), '2c4r')
    cmd.hide('everything')
    cmd.show('cartoon', 'chain L')
    cmd.set('cartoon_transparency', 0.42)
    cmd.color('grey55', 'chain L and polymer')
    cmd.show('sticks', 'chain L and resi 57')
    cmd.show('sticks', 'chain R and resi 10')
    cmd.color('grey40', 'chain L and resi 57')
    cmd.set('stick_radius', 0.26)
    cmd.set_color('rna_blue', RNA_BLUE)
    cmd.color('rna_blue', 'chain R and resi 10 and not name C+*')
    cmd.distance('contact', 'chain L and resi 57 and name CZ', "chain R and resi 10 and name O2'")
    cmd.set('dash_color', 'grey45')
    cmd.set('dash_radius', 0.10)
    cmd.hide('labels', 'contact')
    cmd.set('label_size', 22)
    cmd.set('label_color', 'grey20')
    try:
        cmd.set('label_position', [1.5, 1.5, 1.5])
    except Exception:
        pass
    cmd.orient('resi 57 or resi 10')
    cmd.zoom('resi 57 or resi 10', 4.2)
    cmd.ray(size_w, size_h)
    out = OUT / '2C4R_contact.png'
    cmd.png(str(out), dpi=300)
    print('wrote', out, out.stat().st_size)
    return out


if __name__ == '__main__':
    render('2C4R', PDB / '2C4R.pdb', ['L'], ['R'], zoom_extra=1.5, rna_sticks=True, transparency=0.28)
    render('4QVC', PDB / '4QVC.pdb', ['A', 'B', 'C', 'D', 'E', 'F'], ['G'], chainwise=True, zoom_extra=1.2)
    render('2ZZN', PDB / '2ZZN.pdb', ['A'], ['C'], zoom_extra=1.4, hide_chains=('B', 'D'))
    render('4NGD', ROOT / 'outputs' / 'dicer_input_preparation_20260923' / 'inputs' / 'structures' / '4NGD.pdb',
           ['A'], ['B'], zoom_extra=1.6)
    render_contact()
    print('done ->', OUT)
