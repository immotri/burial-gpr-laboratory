"""Cross-language geometry, heterogeneity, and material slice checks."""
from pathlib import Path
import subprocess,json,sys,numpy as np
root=Path(__file__).resolve().parents[1];sys.path.insert(0,str(root/'toolkit'));import burial_gpr as b
c=b.load(root/'toolkit/model.json')
js="""import fs from 'node:fs';import {geometry,estimate} from './web/core.mjs';import {latentField} from './web/field.mjs';import {sectionRaster} from './web/explorer.mjs';let c=JSON.parse(fs.readFileSync('./toolkit/model.json'));let g=geometry(c),f=latentField(c),slices=[];for(let plane of ['xy','xz','yz']){let extent=plane==='xy'?[.913,3.319,.813,1.783]:plane==='xz'?[.913,3.319,.693,1.337]:[.813,1.783,.693,1.337];let r=sectionRaster(c,g,f,plane,c.sections,extent,161,113);slices.push({plane,extent,data:Array.from(r.data)})}console.log(JSON.stringify({g,f:Array.from(f.data),e:estimate(c),slices}));"""
j=json.loads(subprocess.check_output(['node','--input-type=module','-e',js],cwd=root));g=b.geometry(c);assert len(g)==len(j['g'])
for u,v in zip(g,j['g']):
 for k in u:
  if isinstance(u[k],(float,list)):np.testing.assert_allclose(u[k],v[k],atol=1e-12)
  else:assert u[k]==v[k]
f,h=b.latent(c);np.testing.assert_allclose(f.ravel(),j['f'],atol=1e-12)
for k,v in b.estimate(c).items():np.testing.assert_allclose(v,j['e'][k],atol=1e-12)
for s in j['slices']:
 plane=s['plane'];ext=s['extent'];u=ext[0]+(np.arange(161)+.5)*(ext[1]-ext[0])/161;v=ext[2]+(np.arange(113)+.5)*(ext[3]-ext[2])/113;U,V=np.meshgrid(u,v);p=c['sections']
 X,Y,D=(U,V,np.full_like(U,p['depth_m'])) if plane=='xy' else ((U,np.full_like(U,p['y_m']),V) if plane=='xz' else (np.full_like(U,p['x_m']),U,V))
 ids,_=b.material_at(c,X,Y,D,g,(f,h));np.testing.assert_array_equal(ids.ravel(),s['data'])
print('PASS: detailed geometry, estimates, portable heterogeneity, XY/XZ/YZ material maps agree.')
# Demonstrate that a cranial shell and its sediment-filled cavity are distinct volumes.
c['age_years']=0;g=b.geometry(c);cranium=next(o for o in g if o.get('group')=='cranium' and o['kind']=='ellipsoid');p=np.array(cranium['p']);pts=np.array([p,p+[0,0,-.082]])
a,names=b.material_at(c,pts[:,0],pts[:,1],pts[:,2],g,None);assert [names[x] for x in a]==['fill','bone'],[names[x] for x in a]
print('PASS: cranial cavity and cortical shell are represented separately.')

import tempfile
with tempfile.TemporaryDirectory() as tmp:
 c=b.load(root/'toolkit/smoke_test.json');c['survey'].update(profile_axis='y',polarisation='x',section_offset_m=-.525);c['sections'].update(x_m=1.575,y_m=1.227,depth_m=1.04)
 path=Path(tmp)/'picked.in';meta=b.write_input(c,path)
 assert np.isclose(meta['actual_section_coordinate_m'],1.575)
 assert '#domain: inf ' in path.read_text() and '#hertzian_dipole: x inf ' in path.read_text()
 print('PASS: selected transverse section and correct TM polarisation reach the simulator input.')
