"""General inhumation physics/export regressions, beyond JavaScript/Python parity."""
import copy,tempfile,unittest
from pathlib import Path
import numpy as np
import burial_gpr as b
from cemetery_geometry import contains,bounds,world_point,voxel_boxes,shaft_pieces,frame,grave_state
class CemeteryModels(unittest.TestCase):
 def setUp(self):self.c=b.load(Path(__file__).with_name('cemetery.json'))
 def test_all_starter_scenes_fit(self):
  for path in Path(__file__).with_name('examples').glob('*.json'):
   with self.subTest(path=path.name):b.validate(b.load(path))
 def test_taper_and_step_have_distinct_cross_sections(self):
  g=self.c['graves'][0];s=g['shaft'];s.update(plan='rectangular',inset_m=.15)
  for profile in ['tapered','stepped']:
   s['profile']=profile;objs=[dict(kind='loft',**p) for p in shaft_pieces(s)]
   x=s['length_m']/2-.08;z0=s['top_depth_m']+.03;z1=s['base_depth_m']-.03
   self.assertTrue(any(contains(o,x,0,z0) for o in objs));self.assertFalse(any(contains(o,x,0,z1) for o in objs))
 def test_rotation_changes_world_position_without_shape_change(self):
  g=self.c['graves'][0];g['azimuth_deg']=90;s=shaft_pieces(g['shaft'])[0];o=dict(kind='loft',**s,frame=frame(g));center=[g['x_m'],g['y_m'],.7]
  self.assertTrue(contains(o,*center));self.assertTrue(contains(o,center[0],center[1]+.9,center[2]));self.assertFalse(contains(o,center[0]+.9,center[1],center[2]))
 def test_individual_ages_and_backfill_affect_materials(self):
  g=self.c['graves'][0];h=copy.deepcopy(g);h.update(id='G2',label='G2',seed_offset=992,age_years=0);h['fill']['water_fraction']=.4;g['age_years']=600;self.c['graves'].append(h)
  self.assertAlmostEqual(grave_state(self.c,g)['wood'],.5);self.assertEqual(grave_state(self.c,h)['wood'],1)
  m=b.spectra(self.c,np.array([500e6]));self.assertGreater(m['fill_g1'][0].real,m['fill_g0'][0].real);self.assertNotEqual(m['wood_g0'][0],m['wood_g1'][0])
 def test_recut_order_replaces_earlier_contents(self):
  g=self.c['graves'][0];g['construction']='cist';h=copy.deepcopy(g);h.update(id='G2',label='G2',seed_offset=992,construction='earth');h['contents'].update(bones=False,iron=False,pottery=False);self.c['graves']=[g,h]
  p=[g['x_m'],g['y_m'],g['shaft']['base_depth_m']-.045];a,n=b.material_at(self.c,*p,b.geometry(self.c));self.assertEqual(n[int(a)],'fill_g1')
  self.c['graves'].reverse();a,n=b.material_at(self.c,*p,b.geometry(self.c));self.assertEqual(n[int(a)],'stone')
 def test_non_fitting_container_and_rotated_domain_rejected(self):
  g=self.c['graves'][0];g['construction']='coffin';g['container']['width_m']=2
  with self.assertRaisesRegex(ValueError,'fit'):b.validate(self.c)
  g['construction']='earth';g.update(x_m=.25,azimuth_deg=45)
  with self.assertRaisesRegex(ValueError,'rotated'):b.validate(self.c)
 def test_voxel_box_compression_preserves_exact_occupancy(self):
  c=self.c;g=c['graves'][0];g['azimuth_deg']=27;g['shaft'].update(plan='anthropomorphic',profile='tapered',inset_m=.12);h=.05;ext=(4.8,3.4,2.7);surface=2.4
  # Include a tilted ellipsoid to verify the body frames as well as shaft walls.
  objs=[dict(kind='loft',**shaft_pieces(g['shaft'])[0],frame=frame(g)),dict(kind='ellipsoid',p=[0,0,0],r=[.2,.1,.07],frame={'origin':[2.1,1.5,.8],'rotation':[1,0,0,0,0,-1,0,1,0]})]
  xx=(np.arange(96)+.5)*h;yy=(np.arange(68)+.5)*h;zz=(np.arange(54)+.5)*h;X,Y,Z=np.meshgrid(xx,yy,zz,indexing='ij')
  for o in objs:
   actual=np.zeros(X.shape,bool)
   for a,z in voxel_boxes(o,h,ext,surface):actual|=(X>=a[0])&(X<z[0])&(Y>=a[1])&(Y<z[1])&(Z>=a[2])&(Z<z[2])
   np.testing.assert_array_equal(actual,contains(o,X,Y,surface-Z))
 def test_export_keeps_selected_plane_and_all_materials(self):
  c=self.c;c['survey'].update(cell_m=.03,profile_axis='y',polarisation='x',section_offset_m=.14,dispersive=False);c['soil']['heterogeneity']=0
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'cemetery.in';meta=b.write_input(c,p);self.assertAlmostEqual(meta['actual_section_coordinate_m'],2.54);self.assertIn('fill_g0',p.read_text());self.assertIn('#hertzian_dipole: x inf ',p.read_text())
 def test_preservation_does_not_recentre_surviving_bones(self):
  from cemetery_geometry import body_objects
  g=self.c['graves'][0];a=body_objects(g,dict(bone=1),b.GeometryRNG(42));z=body_objects(g,dict(bone=.4),b.GeometryRNG(42));self.assertTrue(z)
  for o in z:self.assertIn(o,a)
 def test_small_stature_tightens_mesh_screen(self):
  a=b.estimate(self.c)['recommended_cell_m'];self.c['graves'][0]['body']['stature_m']=.7;self.assertLess(b.estimate(self.c)['recommended_cell_m'],a)
