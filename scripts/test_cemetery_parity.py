"""Check all shaft plans/profiles, body poses, grave materials and framed sections."""
from pathlib import Path
import json,subprocess,sys,copy,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'toolkit'))
import burial_gpr as b
from cemetery_geometry import PLANS,PROFILES
models=[]
for path in sorted((ROOT/'toolkit/examples').glob('*.json')):models.append(b.load(path))
base=b.load(ROOT/'toolkit/cemetery.json')
poses=['supine','prone','left_extended','right_extended','left_flexed','right_flexed']
for i,plan in enumerate(PLANS):
 for j,profile in enumerate(PROFILES):
  c=copy.deepcopy(base);g=c['graves'][0];g['shaft'].update(plan=plan,profile=profile,length_m=2.8,width_m=1.5);g['body']['pose']=poses[(i*3+j)%6];g['azimuth_deg']=[13,-31,90][j];c['sections'].update(depth_m=1.14,x_m=2.43,y_m=1.73);c['domain'].update(x_m=5,y_m=4.8);c['soil']['heterogeneity']=0;b.validate(c);models.append(c)
js=json.loads(subprocess.check_output(['node','scripts/cemetery_probe.mjs'],input=json.dumps(models).encode(),cwd=ROOT))
def compare(a,z):
 if isinstance(a,dict):
  assert a.keys()==z.keys(),(a.keys(),z.keys())
  for k in a:compare(a[k],z[k])
 elif isinstance(a,list):
  assert len(a)==len(z)
  for x,y in zip(a,z):compare(x,y)
 elif isinstance(a,(int,float)):np.testing.assert_allclose(a,z,atol=2e-12,rtol=1e-12)
 else:assert a==z,(a,z)
for c,j in zip(models,js):
 compare(b.geometry(c),j['objects'])
 compare(b.estimate(c),j['estimate'])
 for k,v in b.spectra(c,np.array([500e6])).items():np.testing.assert_allclose([v[0].real,v[0].imag],[j['spectra'][k]['r'],j['spectra'][k]['i']],atol=1e-12)
 sig=b.conductivities(c)
 for k,v in j['sigmas'].items():np.testing.assert_allclose(v,sig[k],atol=1e-12)
 # Probe uses no latent field: background heterogeneity remains tested by the legacy script.
 for s in j['slices']:
  ext=s['extent'];u=ext[0]+(np.arange(97)+.5)*(ext[1]-ext[0])/97;v=ext[2]+(np.arange(73)+.5)*(ext[3]-ext[2])/73;U,V=np.meshgrid(u,v);p=c['sections'];plane=s['plane']
  X,Y,D=(U,V,np.full_like(U,p['depth_m'])) if plane=='xy' else ((U,np.full_like(U,p['y_m']),V) if plane=='xz' else (np.full_like(U,p['x_m']),U,V))
  ids,names=b.material_at(c,X,Y,D,b.geometry(c));assert names==s['names'];np.testing.assert_array_equal(ids.ravel(),s['data'])
print(f'PASS: {len(models)} scenes; all five shaft plans, all three wall profiles, all six body positions, materials, estimates and XY/XZ/YZ rasters agree.')
