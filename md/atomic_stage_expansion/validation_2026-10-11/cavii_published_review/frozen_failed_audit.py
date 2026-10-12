"""Compare Ca VII source data with the supplied final Landi-Bhatia article.

Numeric entries were read from published Tables 2-5 and visually checked.
This writes a source audit, not a production atom or an accepted energy model.
"""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
from mocassin_v2_atom_inventory import MOCASSIN,read_atom

# Each row is lower, upper, five-configuration A, twenty-five-configuration A.
TABLE3=[(1,2,.07668,.07638),(1,3,1.620e-5,1.979e-5),(1,4,2.717e-4,1.072e-3),
 (2,3,.1955,.1958),(2,4,1.076,.9256),(2,5,34.74,56.85),(3,4,2.173,1.872),
 (3,5,.1876,.1115),(4,5,4.379,5.186),(3,16,2.623,2.867),(4,16,5.147,5.318),
 (6,16,.4051,.1732),(8,16,1.119,1.396),(9,16,1.056,2.183),(12,16,12.93,11.77),
 (13,16,.01829,.005480),(14,16,5.775e-7,5.905e-7),(15,16,.1961,.1693),
 (2,14,7.806e6,7.075e6),(3,14,4.102e6,3.378e6),(4,14,1.961e7,1.512e7),
 (3,15,1.733e7,1.536e7),(4,15,4.495e6,4.475e6)]
TABLE4=[(1,7,3.038e8,3.092e8),(1,11,8.269e8,4.926e8),(1,18,3.593e9,3.486e9),
 (2,7,1.491e8,1.615e8),(2,8,3.991e8,4.024e8),(2,10,2.518e9,1.508e9),
 (2,11,7.406e8,4.792e8),(2,12,5.327e8,2.928e8),(2,18,1.077e10,9.917e9),
 (3,8,5.396e7,6.468e7),(3,9,4.292e8,4.450e8),(3,11,9.092e8,5.131e8),
 (3,12,1.864e9,1.133e9),(3,18,1.583e10,1.670e10),(4,13,1.838e9,1.825e9),
 (4,17,1.443e10,1.470e10),(5,17,6.632e8,1.602e9)]
TABLE5=[(1,20,1.377e10,1.258e10),(1,23,2.180e10,1.741e10),
 (2,19,2.863e10,2.520e10),(2,20,3.938e9,3.304e9),(2,23,2.277e10,1.918e10),
 (2,24,3.027e10,2.429e10),(3,20,1.076e10,9.271e9),(3,21,1.607e10,1.348e10),
 (3,23,2.602e9,2.331e9),(3,24,1.649e10,1.422e10),(3,25,4.626e10,3.823e10),
 (4,22,4.427e10,3.170e10),(4,26,5.133e10,4.193e10),(4,27,5.262e9,3.909e7),
 (5,27,3.127e10,3.198e10)]
THEORY=[0,1700.388,4265.926,24389.123,50811.305,109858.305,159055.641,
 159127.703,159485.031,185488.125,185531.953,185530.703,204918.625,
 251294.672,252821.609,254954.609,259583.844,250704.281,294410.125,
 293670.312,292080.062,296811.625,301837.781,302536.375,303000.531,
 334259.031,341655.000]

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output-dir',type=Path,required=True);a=p.parse_args()
 out=a.output_dir.resolve()
 if out.exists():raise FileExistsError('Preserve previous published-source review')
 out.mkdir(parents=True)
 article=Path.home()/'MoCHII/references/Landi_2003_ApJ_589_1075.pdf'
 source=MOCASSIN/'data/cavii.dat';atom=read_atom(source,include_tables=True)
 assert atom['nlev']==27
 rates={(r['lower'],r['upper']):r['A_s_inverse'] for r in atom['spontaneous_probabilities']}
 tables={};pairs=[]
 for table,rows,column in [('3',TABLE3,'five_configuration'),('4',TABLE4,'twenty_five_configuration'),('5',TABLE5,'twenty_five_configuration')]:
  records=[]
  for l,u,five,large in rows:
   rate=rates[(l,u)];chosen=five if column=='five_configuration' else large
   err=abs(rate/chosen-1)
   # Printed tables round individual components. Source rows can contain
   # sums of multipoles; agreement here is numeric, not a component identity.
   assert err<.001,(table,l,u,rate,chosen)
   records.append(dict(lower=l,upper=u,source_A_s_inverse=rate,published_five_configuration_A_s_inverse=five,published_twenty_five_configuration_A_s_inverse=large,relative_difference_from_matching_column=err))
   pairs.append((l,u))
  tables[table]=dict(source_matching_column=column,rows=records,maximum_relative_difference=max(r['relative_difference_from_matching_column'] for r in records))
 assert len(set(pairs))==55
 theory=np.array(THEORY);chianti=Path.home()/'RT_Codes/CHIANTI/dbase/ca/ca_7/ca_7.elvlc'
 energies=[]
 for s in chianti.read_text().splitlines():
  if s.strip()=='-1':break
  energies.append(float(s[69:84]))
 assert len(energies)==27
 discrepancy=abs(theory-np.array(energies))
 negative=atom['nonpositive_photon_energy_pairs']
 result=dict(passed=True,final_journal_verified=True,published_pages=10,published_radiation_rows_compared=55,tables=tables,
  published_theoretical_energy_cm_inverse=THEORY,maximum_published_current_theoretical_energy_difference_cm_inverse=float(max(discrepancy)),
  maximum_difference_level=int(np.argmax(discrepancy))+1,source_nonpositive_photon_energy_pairs=negative,
  source_nlev=atom['nlev'],source_positive_radiative_pairs=sum(r['A_s_inverse']>0 for r in atom['spontaneous_probabilities']),
  physical_conclusions=['Printed Table 3 rates agree with the five-configuration column, unlike the adopted twenty-five-configuration radiative model.',
  'Printed Tables 4 and 5 rates agree with the twenty-five-configuration column.',
  'Observed-energy corrections were applied by the authors where both endpoint observations were available.',
  'A theoretical-energy replacement cannot silently retain observed-energy-corrected probabilities as a consistently rescaled model.',
  'Only 55 printed radiation pairs were checked; the configuration provenance of the other source pairs remains unverified.',
  'A total A may include multiple multipoles. Printed M1/E2 values do not always determine the complete total rate.'],
  production_modified=False,corrected_atom_accepted=False,source_sha256={str(f):hashlib.sha256(f.read_bytes()).hexdigest() for f in (Path(__file__),article,source,chianti)})
 (out/'published_source_review.json').write_text(json.dumps(result,indent=2)+'\n')
 print('PASS: final paper and 55 printed radiation rows reviewed; no atom accepted')
if __name__=='__main__':main()
