"""Meaningful numerical and acquisition checks; run: python -m unittest -v."""
import copy,json,tempfile,unittest
from pathlib import Path
import numpy as np
import burial_gpr as b

class PhysicsAndGeometry(unittest.TestCase):
    def setUp(self): self.c=b.load(b.DEFAULT)
    def test_crim_endmembers(self):
        f=np.array([100e6,2e9]);np.testing.assert_allclose(b.crim(f,0,0,4.5),4.5);np.testing.assert_allclose(b.crim(f,1,1,4.5),b.water(f))
    def test_passive_materials_and_fits(self):
        _,rows=b.fit_materials(self.c)
        for name,v in b.spectra(self.c,np.geomspace(50e6,5e9,40)).items():
            self.assertTrue(np.all(v.real>=1),name);self.assertTrue(np.all(v.imag<=1e-12),name)
        for name,v in rows.items():
            self.assertLess(v['max_relative_fit_error'],.01,name);self.assertTrue(all(a>=0 and t>0 for a,t in v['poles']),name)
    def test_decay_limits(self):
        self.c['age_years']=0;self.assertEqual(b.state(self.c)['wood'],1)
        self.c['age_years']=self.c['wood_half_life_years'];self.assertEqual(b.state(self.c)['wood'],.5)
        self.c['environment']='waterlogged';self.assertEqual(b.state(self.c)['void'],0)
    def test_seed_and_frequency_do_not_move_objects(self):
        a=b.geometry(self.c);self.c['survey']['frequency_MHz']=2000;self.assertEqual(a,b.geometry(self.c));self.assertEqual(b.geometry(self.c,True),[])
    def test_invalid_moisture_rejected(self):
        self.c['soil']['water_fraction']=.9
        with self.assertRaises(ValueError):b.validate(self.c)
    def test_geometry_does_not_use_fresh_tissue(self):
        self.assertTrue(all(o['mat'] in ['fill','wood','bone','pottery','pec','rust','air','water'] for o in b.geometry(self.c)))
    def test_control_and_target_share_undisturbed_ground(self):
        field=b.latent(self.c);x=np.array([.1,.2]);y=np.array([.1,.2]);z=np.array([.7,1.5]);a,_=b.material_at(self.c,x,y,z,b.geometry(self.c),field);d,_=b.material_at(self.c,x,y,z,[],field);np.testing.assert_array_equal(a,d)
    def test_geometry_sphere_section(self):
        a,names=b.material_at(self.c,np.array([1.,1.2]),np.array([1.,1.]),np.array([1.,1.]),[dict(kind='sphere',p=[1,1,1],r=.1,mat='bone')]);self.assertEqual(names[a[0]],'bone');self.assertNotEqual(names[a[1]],'bone')
    def test_mesh_guard_checks_whole_sweep(self):
        self.c['survey'].update(cell_m=.0015,frequencies_MHz=[500,2000],frequency_MHz=500)
        with tempfile.TemporaryDirectory() as d:
            with self.assertRaisesRegex(ValueError,'selected sweep'):b.generate(self.c,d)
    def test_acquisition_alignment_and_material_invariance(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'a.in';j=b.write_input(self.c,p);text=p.read_text();self.assertIn('#domain_mode: TM',text)
            q=self.c['survey'];self.assertAlmostEqual(j['trace_spacing_m']/q['cell_m'],round(j['trace_spacing_m']/q['cell_m']))
            self.c['survey']['frequency_MHz']=1000;r,t=b.fit_materials(self.c);self.c['survey']['frequency_MHz']=500;rr,tt=b.fit_materials(self.c);self.assertEqual(r,rr)
    def test_zero_error_for_dense_trace_sampling(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);x=np.arange(11)*.01;t=np.arange(30)*1e-9;a=np.sin(np.arange(30)[:,None]+np.arange(11)[None,:]);np.savez(root/'a.npz',difference=a,x_m=x,time_s=t)
            b.sampling(root/'a.npz',root/'out',[.01,.02],1,25);m=json.loads((root/'out/sampling_metrics.json').read_text());self.assertEqual(m['results'][0]['worst_relative_reconstruction_error'],0);self.assertEqual(m['results'][1]['phases'],2)
    def test_both_grid_sampling_axes_and_phases(self):
        with tempfile.TemporaryDirectory() as d:
            root=Path(d);t=np.arange(8)*1e-9;x=np.arange(9)*.01;y=np.arange(7)*.02;a=np.sin(np.arange(8)[:,None,None]+np.arange(7)[None,:,None]+np.arange(9)[None,None,:]);np.savez(root/'c.npz',difference=a,time_s=t,along_m=x,across_m=y)
            b.grid_sampling(root/'c.npz',root/'out',[.01,.02],[.02,.04],1,6);rows=json.loads((root/'out/grid_metrics.json').read_text())['results'];self.assertEqual(rows[0]['worst_relative_RMS_map_error'],0);self.assertEqual(rows[-1]['phase_combinations'],4)

if __name__=='__main__':unittest.main()
