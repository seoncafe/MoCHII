import json, sys
R = json.load(open(sys.argv[1] if len(sys.argv) > 1 else 'leaf_luminosities.json'))
roman = ['', 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII']
for case, D in R.items():
    rows = D['rows']; tot = D['log']['total']; ch6 = D['log']['ch'][5]
    def new(r): return r.get('rr_fb_exc', 0) + r.get('rr_cas', 0) + r.get('dr_stab', 0) + r.get('dr_cas', 0) + r['L3']
    S = sum(new(r) for r in rows)
    print(f'\n### {case}: total diffuse {tot:.3e} erg/s, channel 6 {ch6:.3e} erg/s\n')
    print('| recombining ion -> product | (1) fb excited | (2) RR cascade | (2) DR stabilization | (2) DR cascade | (3) collisional lines of product... | ground term (Badnell) | new sum |')
    print('|---|---|---|---|---|---|---|---|')
    for r in rows:
        el, st = r['ion'].split('_'); st = int(st)
        n = new(r) - r['L3']
        if n < 0.01 * S: continue
        print(f"| {el.capitalize()}{st-1}+ -> {el.capitalize()} {roman[st-1]} | {r.get('rr_fb_exc',0):.2e} | {r.get('rr_cas',0):.2e} | {r.get('dr_stab',0):.2e} | {r.get('dr_cas',0):.2e} | | {r.get('rr_gt',0):.2e} | {n:.2e} |")
    k = lambda key: sum(r.get(key, 0) for r in rows)
    print(f"| **all** | {k('rr_fb_exc'):.2e} | {k('rr_cas'):.2e} | {k('dr_stab'):.2e} | {k('dr_cas'):.2e} | {k('L3'):.2e} | {k('rr_gt'):.2e} | {S:.2e} |")
    print(f"\nnew/total diffuse = {S/tot:.2e}; new/channel 6 = {S/ch6 if ch6 else float('nan'):.2f}; new/Badnell ground term = {S/k('rr_gt'):.2f}; "
          f"(3)/new = {k('L3')/S:.1e}; He-ionizing (>=24.587 eV) part of new = {(k('rr_fb_exc_he')+k('rr_cas_he')+k('dr_stab_he')+k('dr_cas_he')+k('L3he'))/S:.2f}; "
          f"ionizing share of all metal recombination radiation (ground term + new) / (alpha IP + beta) = {(k('rr_gt')+S-k('L3'))/(k('rr_tot_energy')+k('dr_tot_energy')):.3f}")
    top = sorted(rows, key=lambda r: -r['L3'])[:4]
    print('top (3): ' + ', '.join(f"{r['ion']} {r['L3']:.2e}" for r in top))
