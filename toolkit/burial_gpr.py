#!/usr/bin/env python3
"""Burial GPR Laboratory. Scenario generator, gprMax 4 export and acquisition analysis.
SI units. Age laws are user-controlled hypotheses, not calibrated taphonomy.
"""
from __future__ import annotations
import argparse, copy, csv, hashlib, json, math, os, re, subprocess, sys
from pathlib import Path
import numpy as np
from scipy.ndimage import gaussian_filter
from scipy.optimize import lsq_linear
C = 299792458.0
EPS0 = 8.8541878128e-12
TAUS = np.array([9.231e-12, 5e-11, 2e-10])
VERSION = '1.1.0'
DEFAULT = Path(__file__).with_name('model.json')


def merge(base, change):
    out = copy.deepcopy(base)
    for k,v in change.items():
        out[k] = merge(out[k],v) if isinstance(v,dict) and isinstance(out.get(k),dict) else v
    return out


def load(path):
    c = merge(json.loads(DEFAULT.read_text()), json.loads(Path(path).read_text()))
    validate(c)
    return c


def validate(c):
    def finite(obj):
        if isinstance(obj,dict):
            for v in obj.values(): finite(v)
        elif isinstance(obj,list):
            for v in obj: finite(v)
        elif isinstance(obj,(int,float)) and not math.isfinite(obj): raise ValueError('Non-finite parameter')
    finite(c)
    d,b,s,q=c['domain'],c['coffin'],c['soil'],c['survey']
    for v in [d['x_m'],d['y_m'],d['soil_depth_m'],d['air_m'],b['length_m'],b['width_m'],b['height_m'],b['board_m'],q['cell_m'],q['time_window_ns'],q['trace_spacing_m'],q['line_spacing_m'],s['correlation_m']]:
        if v<=0: raise ValueError('Lengths, spacings and time window must be positive')
    for key,value in c['materials'].items():
        if value<(1 if key.endswith('_eps') else 0): raise ValueError('Material permittivity must be ≥1 and conductivity non-negative')
    if not q['frequencies_MHz']: raise ValueError('Select at least one sweep frequency')
    if c['age_years']<0: raise ValueError('Age must be non-negative')
    for k in ['wood_half_life_years','bone_half_life_years','iron_half_life_years','void_half_life_years']:
        if c[k]<=0: raise ValueError('Half-life hypotheses must be positive')
    if not 0<=c['collapse_limit']<=1: raise ValueError('Collapse limit must be 0–1')
    if c['environment'] not in ['drained','damp','waterlogged']: raise ValueError('Unknown environment')
    for k,w in [('porosity','water_fraction'),('fill_porosity','fill_water_fraction'),('porosity','topsoil_water_fraction')]:
        if not 0<s[k]<1 or not 0<=s[w]<=s[k]: raise ValueError('Water volume fraction must be between 0 and porosity')
    for k in ['wood','bone','pottery']:
        if not 0<s[k+'_porosity']<1 or not 0<=s[k+'_saturation']<=1: raise ValueError('Invalid target porosity/saturation')
    if any(s[k]<0 for k in ['sigma_S_m','fill_sigma_S_m','topsoil_sigma_S_m','heterogeneity']): raise ValueError('Loss and heterogeneity cannot be negative')
    if s['heterogeneity']>0.1: raise ValueError('Water-fraction heterogeneity is capped at 0.1')
    if b['height_m']<.28: raise ValueError('Coffin height must be at least 0.28 m for this object template')
    if b['board_m']>=min(b['width_m'],b['height_m'])/3: raise ValueError('Boards exceed coffin dimensions')
    if b['length_m']<1.85 or b['width_m']<0.58: raise ValueError('This adult skeleton template requires length ≥1.85 m and width ≥0.58 m')
    if b['top_depth_m']<s['topsoil_thickness_m'] or b['top_depth_m']+b['height_m']+0.20>=d['soil_depth_m']: raise ValueError('Coffin must fit beneath topsoil and above the basal margin')
    if b['length_m']+.6>=d['x_m'] or b['width_m']+.6>=d['y_m']: raise ValueError('Add at least 0.3 m domain margin around the coffin')
    if q['mode'] not in ['2D','3D'] or q['profile_axis'] not in ['x','y']: raise ValueError('Invalid simulation mode or profile direction')
    if q['polarisation'] not in ['x','y']: raise ValueError('Use horizontal dipole polarisation x or y')
    if q['mode']=='2D' and q['polarisation']==q['profile_axis']: raise ValueError('2D TM electric dipole must be perpendicular to the profile')
    if q['pml_cells']<8 or int(q['pml_cells'])!=q['pml_cells']: raise ValueError('PML must be an integer of at least 8 cells')
    if q['antenna_height_m']<q['cell_m'] or q['antenna_height_m']+q['cell_m']*(q['pml_cells']+3)>d['air_m']: raise ValueError('Increase air margin or lower antenna height')
    if q['tx_rx_offset_m']<0: raise ValueError('Tx–Rx offset cannot be negative')
    if abs(q['section_offset_m'])>=(d['y_m'] if q['profile_axis']=='x' else d['x_m'])/2-.2: raise ValueError('Section lies outside the usable domain')
    for f in q.get('frequencies_MHz',[])+[q['frequency_MHz']]:
        if not 100<=f<=2000: raise ValueError('Centre frequencies must be 100–2000 MHz')
    for sp in q['spacing_tests_m']:
        if sp<=0: raise ValueError('Sampling tests must be positive')


def state(c):
    age=c['age_years']
    rem=lambda key: 2**(-age/c[key+'_half_life_years'])
    wood=rem('wood')
    return dict(wood=wood,bone=rem('bone'),iron=rem('iron'),void=0 if c['environment']=='waterlogged' else rem('void'),collapse=c['collapse_limit']*(1-wood))


def water(f):
    """20°C, approximate pure-water single Debye pole; excludes DC conductivity."""
    return 4.9+75.2/(1+2j*np.pi*np.asarray(f)*TAUS[0])


def crim(f,porosity,theta,solid):
    return ((1-porosity)*np.sqrt(solid)+theta*np.sqrt(water(f))+(porosity-theta))**2


def spectra(c,f):
    s=c['soil']; st=state(c); m=c['materials']
    vals={
      'subsoil':crim(f,s['porosity'],s['water_fraction'],m['mineral_eps']),
      'topsoil':crim(f,s['porosity'],s['topsoil_water_fraction'],m['mineral_eps']),
      'fill':crim(f,s['fill_porosity'],s['fill_water_fraction'],m['mineral_eps']),
      'bone':crim(f,s['bone_porosity'],s['bone_porosity']*s['bone_saturation'],m['bone_solid_eps']),
      'pottery':crim(f,s['pottery_porosity'],s['pottery_porosity']*s['pottery_saturation'],m['pottery_solid_eps']),
      'water':water(f), 'air':np.ones_like(np.asarray(f),dtype=complex),
      'rust':crim(f,.45,.45*s['wood_saturation'],m['rust_solid_eps'])}
    wood=crim(f,s['wood_porosity'],s['wood_porosity']*s['wood_saturation'],m['wood_solid_eps'])
    vals['wood']=(st['wood']*np.sqrt(wood)+(1-st['wood'])*np.sqrt(vals['fill']))**2
    # Same latent moisture field and same bins at every frequency and in controls.
    for i,t in enumerate(np.linspace(-3,3,13)):
        for base,p,theta in [('subsoil',s['porosity'],s['water_fraction']),('topsoil',s['porosity'],s['topsoil_water_fraction'])]:
            vals[f'{base}_{i}']=crim(f,p,np.clip(theta+t*s['heterogeneity'],0,p),m['mineral_eps'])
    return vals


def conductivities(c):
    s=c['soil']; st=state(c); m=c['materials']
    v=dict(subsoil=s['sigma_S_m'],topsoil=s['topsoil_sigma_S_m'],fill=s['fill_sigma_S_m'],bone=m['bone_sigma_S_m'],pottery=m['pottery_sigma_S_m'],water=m['water_sigma_S_m'],air=0.,rust=m['rust_sigma_S_m'])
    v['wood']=st['wood']*m['wood_sigma_S_m']+(1-st['wood'])*v['fill']
    for i,t in enumerate(np.linspace(-3,3,13)):
        for base in ['subsoil','topsoil']:
            v[f'{base}_{i}']=v[base]*max(.1,1+.15*t)
    return v


def fit_materials(c):
    f=np.geomspace(50e6,5e9,180)
    vals=spectra(c,f); sig=conductivities(c); rows=['#material: 1 0 1 0 air']; tables={}
    basis=np.array([np.ones(f.size),*[1/(1+2j*np.pi*f*t) for t in TAUS]]).T
    A=np.concatenate([basis.real,-basis.imag],axis=0)
    for name,e in vals.items():
        if name=='air': continue
        if c['survey']['dispersive']:
            fit=lsq_linear(A,np.r_[e.real,-e.imag],bounds=([1,0,0,0],[np.inf]*4),tol=1e-12)
            p=fit.x; err=float(np.max(np.abs(basis@p-e)/np.maximum(np.abs(e),1)))
            poles=[(float(a),float(t)) for a,t in zip(p[1:],TAUS) if a>1e-6]
            rows.append(f'#material: {p[0]:.10g} {sig[name]:.10g} 1 0 {name}')
            if poles: rows.append('#add_dispersion_debye: '+str(len(poles))+' '+' '.join(f'{a:.10g} {t:.10g}' for a,t in poles)+' '+name)
            tables[name]=dict(eps_inf=float(p[0]),sigma_dc_S_m=sig[name],poles=poles,max_relative_fit_error=err)
        else:
            e0=complex(spectra(c,np.array([c['survey']['frequency_MHz']*1e6]))[name][0])
            total_sigma=sig[name]+2*np.pi*c['survey']['frequency_MHz']*1e6*EPS0*(-e0.imag)
            rows.append(f'#material: {e0.real:.10g} {total_sigma:.10g} 1 0 {name}')
            tables[name]=dict(eps_r=e0.real,sigma_effective_S_m=total_sigma,approximation='frozen at centre frequency')
    return rows,tables


def estimate(c):
    d,q=c['domain'],c['survey']; h=q['cell_m']; f=q['frequency_MHz']*1e6
    e=spectra(c,np.array([f,2.5*f])); used=['subsoil','topsoil','fill']
    for t in ['wood','bone','pottery']:
        if c['targets'][t if t!='bone' else 'bones']: used.append(t)
    if c['targets']['void'] and state(c)['void']>.15 and c['environment']=='waterlogged': used.append('water')
    # Include wettest moisture bin, because it sets the slowest dielectric wavelength.
    used+=['subsoil_12','topsoil_12']
    if c['targets']['iron']: used.append('rust')
    maxe=max(float(np.max(v.real)) for k,v in e.items() if k in used)
    n=math.sqrt(maxe); wavelength=C/(2.5*f*n)
    extent=d['x_m'] if q['profile_axis']=='x' else d['y_m']
    cells=math.ceil(d['x_m']/h)*math.ceil(d['y_m']/h)*math.ceil((d['soil_depth_m']+d['air_m'])/h) if q['mode']=='3D' else math.ceil(extent/h)*math.ceil((d['soil_depth_m']+d['air_m'])/h)
    dt=.98*h/(C*math.sqrt(3 if q['mode']=='3D' else 2))
    if q['dispersive']: dt=min(dt,TAUS[0]*.45)
    margin=max(.25,(q['pml_cells']+4)*h)
    if 2*margin+q['tx_rx_offset_m']>=extent: raise ValueError('Profile domain is too small for PML and Tx–Rx offset')
    step=max(1,round(q['trace_spacing_m']/h))*h
    ntrace=int(math.floor((extent-2*margin-q['tx_rx_offset_m'])/step))+1
    across=d['y_m'] if q['profile_axis']=='x' else d['x_m']
    line_step=max(1,round(q['line_spacing_m']/h))*h
    nlines=int(math.floor((across-2*margin)/line_step))+1 if q['mode']=='3D' else 1
    detail=min(c['coffin']['board_m'],(.004 if c.get('geometry_version',1)>=2 else .008) if c['targets']['bones'] else .02,.006*math.sqrt(state(c)['iron']) if c['targets']['iron'] else .02,.008 if c['targets']['pottery'] else .02)
    return dict(cells=cells,dt_s=dt,iterations=math.ceil(q['time_window_ns']*1e-9/dt),estimated_RAM_GB_range=[.15+cells*(200 if q['dispersive'] else 80)/1e9,.30+cells*(800 if q['dispersive'] else 200)/1e9],
       wavelength_min_m=wavelength,recommended_cell_m=min(wavelength/15,detail/3),points_per_wavelength=wavelength/h,
       sampling_limit_m=C/(4*2.5*f*math.sqrt(float(e['subsoil'][0].real))),
       ntrace=ntrace,nlines=nlines,actual_trace_spacing_m=step,actual_line_spacing_m=line_step,
       cell_updates_per_trace=cells*math.ceil(q['time_window_ns']*1e-9/dt),
       raw_traces_bytes=ntrace*nlines*math.ceil(q['time_window_ns']*1e-9/dt)*4)


class GeometryRNG:
    # Portable, deterministic xorshift32 for element retention; identical in the UI.
    def __init__(self,seed): self.s=int(seed)&0xffffffff or 1
    def random(self):
        x=self.s; x^=(x<<13)&0xffffffff; x^=x>>17; x^=(x<<5)&0xffffffff; self.s=x&0xffffffff
        return self.s/4294967296.0


def geometry(c,control=False):
    """Ordered boxes, spheres and capped cylinders; physical z is depth-positive.
    Geometry is parametric, not anatomical mesh reconstruction.
    """
    if control: return []
    if c.get('geometry_version',1)>=2:
        from realistic_geometry import geometry_v2
        return geometry_v2(c,state(c),GeometryRNG(c['seed']+201))
    b=c['coffin']; d=c['domain']; st=state(c); tar=c['targets']; rng=GeometryRNG(c['seed']+201)
    xc,yc=d['x_m']/2,d['y_m']/2; L,W,H,t=b['length_m'],b['width_m'],b['height_m'],b['board_m']
    x0,x1=xc-L/2,xc+L/2; y0,y1=yc-W/2,yc+W/2; bottom=b['top_depth_m']+H
    minimum_h=.28 if tar['bones'] else (.24 if tar['pottery'] else 2*t)
    top=bottom-max(H*(1-st['collapse']),minimum_h); floor=bottom-t
    out=[]
    def box(a,b,mat): out.append(dict(kind='box',a=list(a),b=list(b),mat=mat))
    def sphere(p,r,mat): out.append(dict(kind='sphere',p=list(p),r=r,mat=mat))
    def cyl(a,b,r,mat): out.append(dict(kind='cylinder',a=list(a),b=list(b),r=r,mat=mat))
    # Shaft truncates original stratigraphy; greater dimensions than the coffin.
    box((x0-.15,y0-.12,c['soil']['topsoil_thickness_m']),(x1+.15,y1+.12,bottom+.08),'fill')
    box((x0,y0,top),(x1,y1,bottom),'fill')
    # Small surviving cavity is only an explicit hypothesis, and never an ancient body.
    if tar['void'] and st['void']>.01:
        vh=max(.003,(floor-top-t)*st['void'])
        box((x0+t,y0+t,top+t),(x1-t,y1-t,min(floor,top+t+vh)),'water' if c['environment']=='waterlogged' else 'air')
    if tar['wood'] and st['wood']>.002:
        box((x0,y0,bottom-t),(x1,y1,bottom),'wood')
        for a,bb in [((x0,y0,top),(x1,y0+t,bottom)),((x0,y1-t,top),(x1,y1,bottom)),((x0,y0,top),(x0+t,y1,bottom)),((x1-t,y0,top),(x1,y1,bottom))]: box(a,bb,'wood')
        # Depressed lid panels; original geometry remains available through age=0.
        for i in range(12):
            xa=x0+L*i/12; xb=x0+L*(i+1)/12
            drop=min((floor-top-t)*.5,.035)*st['collapse']*math.sin(math.pi*(i+.5)/12)**2
            box((xa,y0,top+drop),(xb,y1,top+t+drop),'wood')
    # Scale longitudinal locations to fit editable coffin length.
    xx=lambda u: x0+.12+u*(L-.24)/1.85
    bz=min(floor-.04,max(top+.12,floor-.08))
    if tar['bones']:
        # Seeded retention is correlated by whole anatomical element, not voxels.
        def bonec(a,b,r):
            if rng.random()<st['bone']: cyl(a,b,r,'bone')
        if rng.random()<st['bone']:
            sphere((xx(1.66),yc,bz-.035),.085,'bone'); sphere((xx(1.66),yc,bz-.035),.077,'fill')
        bonec((xx(.86),yc,bz),(xx(1.48),yc,bz),.017)
        for side in [-1,1]:
            yy=yc+side*.095
            bonec((xx(.13),yc+side*.075,bz),(xx(.53),yc+side*.10,bz),.013)
            bonec((xx(.53),yc+side*.10,bz),(xx(.94),yc+side*.09,bz),.017)
            bonec((xx(.97),yc+side*.05,bz),(xx(1.02),yc+side*.12,bz),.035)
            bonec((xx(1.02),yc+side*.16,bz),(xx(1.23),yc+side*.19,bz),.011)
            bonec((xx(1.23),yc+side*.19,bz),(xx(1.47),yc+side*.16,bz),.013)
            for u in np.linspace(1.12,1.43,7):
                # Rib half-arcs made of short cylinders, 8 mm diameter.
                if rng.random()>st['bone']: continue
                pts=[(xx(u),yc+side*.13*math.sin(a),bz-.055*math.sin(a)) for a in np.linspace(0,math.pi/2,5)]
                for a,bb in zip(pts[:-1],pts[1:]): cyl(a,bb,.004,'bone')
    if tar['iron']:
        # Corrosion halo and a shrinking PEC core: high-conductivity limit, no ferromagnetic law.
        sx,sy=xx(1.04),yc-.245
        box((sx-.41,sy-.025,bz-.018),(sx+.41,sy+.025,bz+.018),'rust')
        if st['iron']>.005:
            thick=.006*math.sqrt(st['iron']); width=.035*math.sqrt(st['iron'])
            box((sx-.40,sy-width/2,bz-thick/2),(sx+.40,sy+width/2,bz+thick/2),'pec')
            box((sx+.27,sy-.07,bz-thick/2),(sx+.30,sy+.07,bz+thick/2),'pec')
    if tar['pottery']:
        for px,py,r in [(xx(.13),yc-.14,.065),(xx(.23),yc+.13,.07),(xx(.05),yc+.02,.052)]:
            pz=floor-r
            sphere((px,py,pz),r,'pottery'); sphere((px,py,pz),r-.008,'fill')
            # Open neck: vessel interior is sediment-filled by default.
            cyl((px,py,pz-r-.01),(px,py,pz-r+.035),r*.35,'fill')
    return out


def latent(c):
    d=c['domain']; corr=c['soil']['correlation_m']; step=corr/2
    shape=tuple(max(3,int(math.ceil(d[k]/step))+1) for k in ['x_m','y_m','soil_depth_m'])
    if np.prod(shape)>15_000_000: raise ValueError('Correlation too short for this domain')
    if c.get('geometry_version',1)>=2:
        rng=GeometryRNG(c['seed']+8000); modes=[]
        for _ in range(12): modes.append([*[(rng.random()*2-1)*math.pi/corr for _ in range(3)],rng.random()*2*math.pi])
        I,J,K=np.indices(shape,dtype=float);a=np.zeros(shape)
        for x,y,z,phase in modes:a+=np.sin(x*I*step+y*J*step+z*K*step+phase)
    else: a=gaussian_filter(np.random.default_rng(c['seed']).normal(size=shape),1.0,mode='reflect')
    a=(a-a.mean())/max(a.std(),1e-9)
    return a,step


def material_at(c,X,Y,D,objects,field=None):
    X,Y,D=np.broadcast_arrays(X,Y,D); names=['air','subsoil','topsoil','fill','wood','bone','pottery','water','rust','pec']+[f'{base}_{i}' for base in ['subsoil','topsoil'] for i in range(13)]
    ids={k:i for i,k in enumerate(names)}; a=np.full(X.shape,ids['subsoil'],np.int16)
    a[D<0]=ids['air']; a[(D>=0)&(D<c['soil']['topsoil_thickness_m'])]=ids['topsoil']
    if field is not None and c['soil']['heterogeneity']>0:
        v,step=field; ix=np.clip((X/step).astype(int),0,v.shape[0]-1); iy=np.clip((Y/step).astype(int),0,v.shape[1]-1); iz=np.clip((D/step).astype(int),0,v.shape[2]-1)
        bins=np.rint(np.clip(v[ix,iy,iz],-3,3)*2+6).astype(int)
        for base in ['subsoil','topsoil']:
            mask=a==ids[base]; a[mask]=ids[base+'_0']+bins[mask]
    bounds=[(float(v.min()),float(v.max())) for v in [X,Y,D]]
    for o in objects:
        if o['kind']=='box':lo,hi=o['a'],o['b']
        elif o['kind'] in ['sphere','ellipsoid']:
            rr=[o['r']]*3 if o['kind']=='sphere' else o['r'];lo=[o['p'][k]-rr[k] for k in range(3)];hi=[o['p'][k]+rr[k] for k in range(3)]
        else:lo=[min(o['a'][k],o['b'][k])-o['r'] for k in range(3)];hi=[max(o['a'][k],o['b'][k])+o['r'] for k in range(3)]
        if any(bounds[k][1]<lo[k] or bounds[k][0]>hi[k] for k in range(3)):continue
        if o['kind']=='box':
            m=(X>=o['a'][0])&(X<o['b'][0])&(Y>=o['a'][1])&(Y<o['b'][1])&(D>=o['a'][2])&(D<o['b'][2])
        elif o['kind']=='ellipsoid':m=((X-o['p'][0])/o['r'][0])**2+((Y-o['p'][1])/o['r'][1])**2+((D-o['p'][2])/o['r'][2])**2<=1
        elif o['kind']=='sphere': m=(X-o['p'][0])**2+(Y-o['p'][1])**2+(D-o['p'][2])**2<=o['r']**2
        else:
            p=np.array(o['a']); v=np.array(o['b'])-p; vv=float(v@v)
            u=((X-p[0])*v[0]+(Y-p[1])*v[1]+(D-p[2])*v[2])/vv
            m=(u>=0)&(u<=1)&((X-p[0]-u*v[0])**2+(Y-p[1]-u*v[1])**2+(D-p[2]-u*v[2])**2<=o['r']**2)
        a[m]=ids[o['mat']]
    return a,names


def snap(v,h): return round(v/h)*h

def write_input(c,path,control=False):
    d,q=c['domain'],c['survey']; h=q['cell_m']; est=estimate(c); rows,tables=fit_materials(c)
    # Explicit snap avoids different domain rounding in solver and generator.
    nx,ny,nz=[math.ceil(v/h) for v in [d['x_m'],d['y_m'],d['soil_depth_m']+d['air_m']]]
    ext=(nx*h,ny*h,nz*h); surface=snap(d['soil_depth_m'],h)
    mode=q['mode']; along=q['profile_axis']; invariant='y' if along=='x' else 'x'
    if mode=='2D': domain=f'{ext[0]:.10g} inf {ext[2]:.10g}' if along=='x' else f'inf {ext[1]:.10g} {ext[2]:.10g}'
    else: domain=' '.join(f'{v:.10g}' for v in ext)
    dtfull=h/(C*math.sqrt(2 if mode=='2D' else 3)); factor=min(.98,TAUS[0]*.45/dtfull) if q['dispersive'] else .98
    lines=['#title: Burial GPR Laboratory '+('background control' if control else 'burial target'),f'#domain: {domain}',f'#dx_dy_dz: {h:.10g} {h:.10g} {h:.10g}',f'#time_window: {q["time_window_ns"]*1e-9:.10e}',f'#time_step_stability_factor: {factor:.10g}']
    if mode=='2D': lines.append('#domain_mode: TM')
    p=q['pml_cells']; lines.append(f'#pml_cells: {p if invariant!="x" or mode=="3D" else 0} {p if invariant!="y" or mode=="3D" else 0} {p} {p if invariant!="x" or mode=="3D" else 0} {p if invariant!="y" or mode=="3D" else 0} {p}')
    lines+=rows
    omitted=[]
    objs=geometry(c,control); field=latent(c) if c['soil']['heterogeneity']>0 else None
    def box(a,b,mat):
        # All boxes generated on explicit voxel edges; no subcell expansion.
        if all(bb>aa for aa,bb in zip(a,b)):
            lines.append('#box: '+' '.join(f'{v:.10g}' for v in [*a,*b])+' '+mat+' n')
    if mode=='2D':
        n=nx if along=='x' else ny; u=(np.arange(n)+.5)*h; z=(np.arange(nz)+.5)*h; dep=surface-z
        U,D=np.meshgrid(u,dep,indexing='ij')
        section=(d['y_m']/2 if along=='x' else d['x_m']/2)+q['section_offset_m']
        X,Y=(U,np.full_like(U,section)) if along=='x' else (np.full_like(U,section),U)
        ids,names=material_at(c,X,Y,D,objs,field)
        for j in range(nz):
            col=ids[:,j]; breaks=np.r_[0,np.flatnonzero(np.diff(col))+1,n]
            for a,b in zip(breaks[:-1],breaks[1:]):
                mat=names[col[a]]
                if mat=='air': continue
                if along=='x': box((a*h,0,j*h),(b*h,h,(j+1)*h),mat)
                else: box((0,a*h,j*h),(h,b*h,(j+1)*h),mat)
    else:
        box((0,0,0),(ext[0],ext[1],surface),'subsoil'); top=snap(c['soil']['topsoil_thickness_m'],h)
        box((0,0,surface-top),(ext[0],ext[1],surface),'topsoil')
        if field is not None:
            fld,step=field
            # Tiles are coarser than the electromagnetic mesh and fixed in physical space.
            for i in range(math.ceil(ext[0]/step)):
                for j in range(math.ceil(ext[1]/step)):
                    for k in range(math.ceil(surface/step)):
                        dep0=k*step; dep1=min(surface,(k+1)*step)
                        val=fld[min(i,fld.shape[0]-1),min(j,fld.shape[1]-1),min(k,fld.shape[2]-1)]
                        bi=int(round(np.clip(val,-3,3)*2+6))
                        # Split tile at the stratigraphic boundary.
                        cuts=sorted(set([dep0,dep1]+([c['soil']['topsoil_thickness_m']] if dep0<c['soil']['topsoil_thickness_m']<dep1 else [])))
                        for da,db in zip(cuts[:-1],cuts[1:]):
                            base='topsoil' if (da+db)/2<c['soil']['topsoil_thickness_m'] else 'subsoil'
                            a=[snap(i*step,h),snap(j*step,h),snap(surface-db,h)]; b=[min(ext[0],snap((i+1)*step,h)),min(ext[1],snap((j+1)*step,h)),snap(surface-da,h)]
                            box(a,b,base+'_'+str(bi))
        def coord(p): return [snap(p[0],h),snap(p[1],h),snap(surface-p[2],h)]
        for o in objs:
            if o['kind']=='box':
                a=coord(o['a']); b=coord(o['b']); lo=[min(x,y) for x,y in zip(a,b)]; hi=[max(x,y) for x,y in zip(a,b)]; box(lo,hi,o['mat'])
            elif o['kind']=='ellipsoid':
                if min(o['r'])<h:
                    omitted.append(dict(reason='ellipsoid minor semi-axis smaller than a cell',object=o));continue
                lines.append('#ellipsoid: '+' '.join(f'{v:.10g}' for v in [*coord(o['p']),*o['r']])+' '+o['mat']+' n')
            elif o['kind']=='sphere':
                if o['r']<h:
                    omitted.append(dict(reason='radius smaller than one mesh cell',object=o)); continue
                lines.append('#sphere: '+' '.join(f'{v:.10g}' for v in [*coord(o['p']),o['r']])+' '+o['mat']+' n')
            else:
                a,b=coord(o['a']),coord(o['b'])
                if o['r']<h or a==b:
                    omitted.append(dict(reason='subcell radius or coincident snapped endpoints',object=o)); continue
                lines.append('#cylinder: '+' '.join(f'{v:.10g}' for v in [*a,*b,o['r']])+' '+o['mat']+' n')
    margin=max(.25,(p+4)*h); start=snap(margin,h); zsrc=surface+snap(q['antenna_height_m'],h); offset=snap(q['tx_rx_offset_m'],h)
    across=ext[1] if along=='x' else ext[0]
    cross=snap(across/2+q['section_offset_m'],h)
    source=[start,cross,zsrc] if along=='x' else [cross,start,zsrc]
    recv=source.copy(); recv[0 if along=='x' else 1]+=offset
    if mode=='2D': source[1 if along=='x' else 0]='inf'; recv[1 if along=='x' else 0]='inf'
    point=lambda a: ' '.join(v if isinstance(v,str) else f'{v:.10g}' for v in a)
    step=[est['actual_trace_spacing_m'],0,0] if along=='x' else [0,est['actual_trace_spacing_m'],0]
    lines += [f'#waveform: ricker 1 {q["frequency_MHz"]*1e6:.10g} pulse',f'#hertzian_dipole: {q["polarisation"]} {point(source)} pulse',f'#rx: {point(recv)}',f'#src_steps: {point(step)}',f'#rx_steps: {point(step)}']
    path.write_text('\n'.join(lines)+'\n')
    return dict(file=path.name,ntraces=est['ntrace'],actual_section_coordinate_m=(section if mode=='2D' else cross),frequency_MHz=q['frequency_MHz'],mode=mode,profile_axis=along,section_offset_m=q['section_offset_m'],start_m=start,trace_spacing_m=est['actual_trace_spacing_m'],tx_rx_offset_m=offset,receiver_component='E'+q['polarisation'],materials=tables,estimate=est,control=control,omitted_subcell_primitives=omitted)


def plot_geometry(c,out):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap
    objs=geometry(c);field=latent(c);d=c['domain'];p=c.get('sections',dict(x_m=d['x_m']/2,y_m=d['y_m']/2,depth_m=c['coffin']['top_depth_m']+c['coffin']['height_m']-.09))
    fig,axes=plt.subplots(1,3,figsize=(15,5),layout='constrained')
    colors=['#dce9e4','#bcab8d','#737856','#c9ad88','#765038','#e6dbc0','#aa6041','#669ba5','#966b42','#3c4140']
    def tint(hexcol,factor):return '#'+''.join(f'{min(255,max(0,round(int(hexcol[i:i+2],16)*factor))):02x}' for i in [1,3,5])
    colors += [tint(base,1-(i-6)*.025) for base in ['#bcab8d','#737856'] for i in range(13)]
    cmap=ListedColormap(colors)
    for ax,plane in zip(axes,['xy','xz','yz']):
        u=np.linspace(0,d['y_m'] if plane=='yz' else d['x_m'],650);v=np.linspace(0,d['y_m'] if plane=='xy' else d['soil_depth_m'],420);U,V=np.meshgrid(u,v)
        X,Y,D=(U,V,np.full_like(U,p['depth_m'])) if plane=='xy' else ((U,np.full_like(U,p['y_m']),V) if plane=='xz' else (np.full_like(U,p['x_m']),U,V))
        ids,names=material_at(c,X,Y,D,objs,field)
        ax.imshow(ids,origin='upper',extent=[u[0],u[-1],v[-1],v[0]],vmin=0,vmax=len(names)-1,cmap=cmap,interpolation='nearest')
        title=f'XY at depth {p["depth_m"]:.3f} m' if plane=='xy' else (f'XZ at y={p["y_m"]:.3f} m' if plane=='xz' else f'YZ at x={p["x_m"]:.3f} m')
        ax.set(title=title,xlabel='y (m)' if plane=='yz' else 'x (m)',ylabel='y (m)' if plane=='xy' else 'Depth (m)')
        ax.axvline(p['y_m'] if plane=='yz' else p['x_m'],color='#477c99',lw=.7,ls='--');ax.axhline(p['y_m'] if plane=='xy' else p['depth_m'],color='#b58b2b',lw=.7,ls='--')
    fig.suptitle(f'Geometry v{c.get("geometry_version",1)} | {c["age_years"]} years | {c["environment"]} hypothesis | material sections, not radar data');fig.savefig(out,dpi=160);plt.close(fig)


def generate(c,out,allow=False,raster=False):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    est=estimate(c)
    strict=copy.deepcopy(c); strict['survey']['frequency_MHz']=max(c['survey']['frequencies_MHz'])
    worst=estimate(strict)
    if c['survey']['cell_m']>worst['recommended_cell_m']*1.02 and not allow:
        raise ValueError(f'Mesh is under-resolved: requested {c["survey"]["cell_m"]*1000:.2f} mm, recommendation ≤{worst["recommended_cell_m"]*1000:.2f} mm for the selected sweep. Refine or use --allow-underresolved for a labelled screening experiment.')
    ctext=json.dumps(c,indent=2)
    digest=hashlib.sha256(ctext.encode()).hexdigest()
    if (out/'manifest.json').exists() and json.loads((out/'manifest.json').read_text()).get('project_sha256')!=digest: raise ValueError('This output folder belongs to a different project. Choose a new --out folder to keep results separate.')
    (out/'project.json').write_text(ctext)
    manifest=dict(tool_version=VERSION,generator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),project_sha256=digest,scientific_status='scenario hypotheses; ideal dipole; not calibrated to a commercial antenna',underresolved=c['survey']['cell_m']>worst['recommended_cell_m']*1.02,jobs=[])
    for f in c['survey']['frequencies_MHz']:
        cf=copy.deepcopy(c); cf['survey']['frequency_MHz']=f
        offsets=[c['survey']['section_offset_m']]
        if raster and c['survey']['mode']=='3D':
            across=c['domain']['y_m' if c['survey']['profile_axis']=='x' else 'x_m']; m=max(.25,(c['survey']['pml_cells']+4)*c['survey']['cell_m']); step=estimate(cf)['actual_line_spacing_m']
            offsets=(np.arange(estimate(cf)['nlines'])*step+m-across/2).tolist()
        for li,off in enumerate(offsets):
            cf['survey']['section_offset_m']=off
            for control in [False,True]:
                name=f'f{f:g}_line{li:03d}_'+('background' if control else 'burial')+'.in'
                job=write_input(cf,out/name,control); job['underresolved']=cf['survey']['cell_m']>job['estimate']['recommended_cell_m']*1.02; manifest['jobs'].append(job)
    (out/'manifest.json').write_text(json.dumps(manifest,indent=2)); plot_geometry(c,out/'geometry.png')
    print(json.dumps(dict(directory=str(out.resolve()),jobs=len(manifest['jobs']),estimate=est,underresolved=manifest['underresolved']),indent=2))


def run_jobs(folder,limit,threads,first_job=False):
    root=Path(folder).resolve(); man=json.loads((root/'manifest.json').read_text())
    for job in man['jobs'][:1] if first_job else man['jobs']:
        n=min(job['ntraces'],limit) if limit else job['ntraces']
        cmd=[sys.executable,'-m','gprMax',job['file'],'-n',str(n),'--geometry-fixed']
        if job['underresolved']: cmd.append('--allow-underresolved')
        env=dict(os.environ); env['OMP_NUM_THREADS']=str(threads)
        print('Running',job['file'],n,'traces',flush=True)
        subprocess.run(cmd,cwd=root,env=env,check=True)
    print('Finished requested jobs. Collect only complete background/target pairs for analysis.')


def read_scan(root,job):
    import h5py
    stem=Path(job['file']).stem; files=[]
    for p in root.glob(stem+'*.h5'):
        suffix=p.stem[len(stem):]
        if suffix=='' or suffix.isdigit(): files.append((int(suffix or 1),p))
    files=sorted(files)
    if not files: raise FileNotFoundError('No numbered gprMax .h5 outputs for '+stem)
    traces=[]; dt0=None; indices=[]
    for ix,path in files:
        with h5py.File(path,'r') as f:
            dt=float(f.attrs['dt']); keys=list(f['rxs'].keys())
            if len(keys)!=1: raise ValueError('Expected one receiver per trace; verify identity')
            rx=f['rxs'][keys[0]]; data=np.asarray(rx[job['receiver_component']])
            if 'Position' in rx.attrs:
                actual=np.array(rx.attrs['Position']); along=0 if job['profile_axis']=='x' else 1
                expected=job['start_m']+(ix-1)*job['trace_spacing_m']+job['tx_rx_offset_m']
                if not np.isclose(actual[along],expected,atol=1e-7): raise ValueError('Receiver position does not match acquisition manifest')
            if data.ndim!=1: raise ValueError('Expected unmerged one-dimensional receiver histories')
            if dt0 is not None and not np.isclose(dt,dt0,rtol=1e-9,atol=0): raise ValueError('Trace time steps differ')
            traces.append(data); indices.append(ix); dt0=dt
    if indices!=list(range(1,len(files)+1)): raise ValueError('Trace sequence has missing indices')
    return np.stack(traces,axis=1),dt0


def collect(folder):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root=Path(folder); man=json.loads((root/'manifest.json').read_text()); jobs=man['jobs']; pairs={}
    for job in jobs: pairs.setdefault((job['frequency_MHz'],job['section_offset_m']),{})['background' if job['control'] else 'burial']=job
    produced={}
    for (freq,section),p in pairs.items():
        if len(p)!=2: continue
        try: a,dt=read_scan(root,p['burial']); b,bdt=read_scan(root,p['background'])
        except FileNotFoundError as exc: print('Skipping incomplete pair:',exc); continue
        if a.shape!=b.shape or not np.isclose(dt,bdt,rtol=1e-9,atol=0): raise ValueError('Background and burial acquisitions do not match')
        job=p['burial']; x=job['start_m']+np.arange(a.shape[1])*job['trace_spacing_m']; t=np.arange(a.shape[0])*dt
        diff=a-b; stem=Path(job['file']).stem
        np.savez_compressed(root/(stem+'_scan.npz'),raw=a,background=b,difference=diff,time_s=t,x_m=x,metadata=json.dumps(job))
        if job['mode']=='3D': produced.setdefault(freq,[]).append((job.get('actual_section_coordinate_m',section),root/(stem+'_scan.npz')))

        fig,axs=plt.subplots(1,3,figsize=(13,5),layout='constrained'); clip=np.max(np.abs(np.r_[a.ravel(),b.ravel(),diff.ravel()])) or 1
        for ax,data,title in zip(axs,[a,b,diff],['Burial','Background','Burial − background']):
            ax.imshow(data,aspect='auto',cmap='RdBu',vmin=-clip,vmax=clip,extent=[x[0]-job['trace_spacing_m']/2,x[-1]+job['trace_spacing_m']/2,(t[-1]+dt/2)*1e9,-dt/2*1e9]); ax.set(xlabel='Tx position (m)',ylabel='Time (ns)',title=title)
        fig.suptitle(f'gprMax full-wave {job["mode"]}, {freq:g} MHz | same amplitude scale | '+('UNDER-RESOLVED SCREENING' if job['underresolved'] else 'convergence still required'))
        fig.savefig(root/(stem+'_radar.png'),dpi=160); plt.close(fig)
    for freq,lines in produced.items():
        lines.sort(); expected=sum(1 for job in jobs if not job['control'] and job['frequency_MHz']==freq)
        if len(lines)<2 or len(lines)!=expected: continue
        arrays=[]; template=None
        for cross,file in lines:
            scan=np.load(file)
            if template is not None and (scan['difference'].shape!=template['difference'].shape or not np.allclose(scan['x_m'],template['x_m']) or not np.allclose(scan['time_s'],template['time_s'])): raise ValueError('Raster lines have incompatible timing or sampling')
            arrays.append(scan['difference']); template=scan
        np.savez_compressed(root/f'f{freq:g}_raster.npz',difference=np.stack(arrays,axis=1),time_s=template['time_s'],along_m=template['x_m'],across_m=np.array([cross for cross,file in lines]))
    print('Saved NPZ arrays and radar PNGs for available matched pairs; complete 3D rasters also saved as time × line × trace cubes.')


def sampling(scan,out,spacings,gate_start_ns=0,gate_end_ns=None):
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    f=np.load(scan); data=f['difference']; x=f['x_m']; t=f['time_s']; base=float(np.median(np.diff(x)))
    if len(x)<3: raise ValueError('At least three traces are needed for sampling analysis')
    gate=t>=gate_start_ns*1e-9
    if gate_end_ns is not None: gate &= t<=gate_end_ns*1e-9
    if not gate.any(): raise ValueError('Empty analysis time gate')
    target=data[gate]; denom=np.linalg.norm(target)
    fig,axes=plt.subplots(1,len(spacings),figsize=(4*len(spacings),5),layout='constrained'); axes=np.atleast_1d(axes)
    rows=[]; clip=np.max(np.abs(data)) or 1
    for ax,spacing in zip(axes,spacings):
        stride=round(spacing/base)
        if stride<1 or not np.isclose(stride*base,spacing,atol=1e-7): raise ValueError('Each spacing must be an integer multiple of actual dense spacing')
        errors=[]
        for phase in range(min(stride,len(x)-1)):
            ix=np.arange(phase,len(x),stride)
            if len(ix)<2: continue
            # Only compare within sampled support; do not extrapolate endpoints.
            dense=np.arange(ix[0],ix[-1]+1); recon=np.stack([np.interp(x[dense],x[ix],row[ix]) for row in target])
            den=np.linalg.norm(target[:,dense]); errors.append(float(np.linalg.norm(recon-target[:,dense])/den) if den else float('nan'))
        ix=np.arange(0,len(x),stride)
        ax.imshow(data[:,ix],aspect='auto',interpolation='nearest',cmap='RdBu',vmin=-clip,vmax=clip,extent=[x[ix[0]]-spacing/2,x[ix[-1]]+spacing/2,t[-1]*1e9,0]); ax.set(title=f'{spacing*100:.1f} cm | {len(ix)} traces',xlabel='Tx position (m)',ylabel='Time (ns)',xlim=(x[0]-base/2,x[-1]+base/2))
        rows.append(dict(spacing_m=spacing,stride=stride,phases=len(errors),best_relative_reconstruction_error=min(errors) if errors else None,worst_relative_reconstruction_error=max(errors) if errors else None))
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    meta=json.loads(str(f['metadata'])) if 'metadata' in f else {}
    warning='\nUNDER-RESOLVED SCREENING: numerical accuracy is not established' if meta.get('underresolved') else ''
    fig.suptitle('Unfiltered trace thinning | phase 0 shown | common amplitude scale'+warning); fig.savefig(out/'sampling.png',dpi=160); plt.close(fig)
    (out/'sampling_metrics.json').write_text(json.dumps(dict(time_gate_ns=[gate_start_ns,gate_end_ns],metric='relative L2 error after linear spatial interpolation, same sampled support; not detection probability',results=rows),indent=2))
    print(json.dumps(rows,indent=2))


def grid_sampling(cube,out,trace_spacings,line_spacings,gate_start_ns,gate_end_ns):
    from scipy.interpolate import RegularGridInterpolator
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    f=np.load(cube); data=f['difference']; t=f['time_s']; x=f['along_m']; y=f['across_m']
    if data.ndim!=3 or len(x)<3 or len(y)<3: raise ValueError('Need a complete 3D raster with ≥3 traces and ≥3 lines')
    gate=(t>=gate_start_ns*1e-9)&(t<=gate_end_ns*1e-9)
    if not gate.any(): raise ValueError('Empty time gate')
    # RMS time-slice is a reproducible map statistic, not a migrated depth slice.
    amplitude=np.sqrt(np.mean(data[gate]**2,axis=0)); dx=float(np.median(np.diff(x))); dy=float(np.median(np.diff(y)))
    results=[]; fig,axes=plt.subplots(len(line_spacings),len(trace_spacings),squeeze=False,figsize=(4*len(trace_spacings),3.5*len(line_spacings)),layout='constrained')
    for i,ls in enumerate(line_spacings):
        ky=round(ls/dy)
        if ky<1 or not np.isclose(ky*dy,ls,atol=1e-7): raise ValueError('Line spacings must be integer multiples of actual dense line spacing')
        for j,ts in enumerate(trace_spacings):
            kx=round(ts/dx)
            if kx<1 or not np.isclose(kx*dx,ts,atol=1e-7): raise ValueError('Trace spacings must be integer multiples of actual dense trace spacing')
            errors=[]
            for py in range(min(ky,len(y)-1)):
                iy=np.arange(py,len(y),ky)
                if len(iy)<2: continue
                for px in range(min(kx,len(x)-1)):
                    ix=np.arange(px,len(x),kx)
                    if len(ix)<2: continue
                    xx=np.arange(ix[0],ix[-1]+1); yy=np.arange(iy[0],iy[-1]+1); Y,X=np.meshgrid(y[yy],x[xx],indexing='ij')
                    sampled=amplitude[np.ix_(iy,ix)]; recon=RegularGridInterpolator((y[iy],x[ix]),sampled)(np.stack([Y,X],axis=-1)); truth=amplitude[np.ix_(yy,xx)]; den=np.linalg.norm(truth)
                    errors.append(float(np.linalg.norm(recon-truth)/den) if den else float('nan'))
            ix=np.arange(0,len(x),kx); iy=np.arange(0,len(y),ky); ax=axes[i,j]
            ax.imshow(amplitude[np.ix_(iy,ix)],origin='lower',aspect='equal',interpolation='nearest',vmin=0,vmax=float(amplitude.max()) or 1,cmap='magma',extent=[x[ix[0]]-ts/2,x[ix[-1]]+ts/2,y[iy[0]]-ls/2,y[iy[-1]]+ls/2]); ax.set(title=f'Trace {ts*100:g} cm / line {ls*100:g} cm',xlabel='Along profile (m)',ylabel='Across profiles (m)',xlim=(x[0]-dx/2,x[-1]+dx/2),ylim=(y[0]-dy/2,y[-1]+dy/2))
            results.append(dict(trace_spacing_m=ts,line_spacing_m=ls,phase_combinations=len(errors),best_relative_RMS_map_error=min(errors) if errors else None,worst_relative_RMS_map_error=max(errors) if errors else None))
    out=Path(out);out.mkdir(parents=True,exist_ok=True);fig.suptitle(f'RMS difference amplitude, {gate_start_ns:g}–{gate_end_ns:g} ns | no migration');fig.savefig(out/'grid_sampling.png',dpi=160);plt.close(fig)
    (out/'grid_metrics.json').write_text(json.dumps(dict(metric='relative L2 error of linearly reconstructed gated RMS maps, all starting phases; not detection probability',results=results),indent=2));print(json.dumps(results,indent=2))


def main():
    parser=argparse.ArgumentParser(description=__doc__); sub=parser.add_subparsers(dest='cmd',required=True)
    p=sub.add_parser('estimate'); p.add_argument('project')
    p=sub.add_parser('generate'); p.add_argument('project'); p.add_argument('--out',default='runs'); p.add_argument('--allow-underresolved',action='store_true'); p.add_argument('--raster',action='store_true')
    p=sub.add_parser('run'); p.add_argument('folder'); p.add_argument('--limit-traces',type=int); p.add_argument('--threads',type=int,default=4); p.add_argument('--first-job',action='store_true')
    p=sub.add_parser('collect'); p.add_argument('folder')
    p=sub.add_parser('sampling'); p.add_argument('scan'); p.add_argument('--out',default='sampling'); p.add_argument('--spacings',nargs='+',type=float,default=[.01,.02,.05,.10]); p.add_argument('--gate-start-ns',type=float,default=0); p.add_argument('--gate-end-ns',type=float)
    p=sub.add_parser('grid-sampling'); p.add_argument('cube'); p.add_argument('--out',default='grid-sampling'); p.add_argument('--trace-spacings',nargs='+',type=float,required=True); p.add_argument('--line-spacings',nargs='+',type=float,required=True); p.add_argument('--gate-start-ns',type=float,required=True); p.add_argument('--gate-end-ns',type=float,required=True)
    args=parser.parse_args()
    try:
        if args.cmd=='estimate': print(json.dumps(dict(state=state(load(args.project)),estimate=estimate(load(args.project))),indent=2))
        elif args.cmd=='generate': generate(load(args.project),args.out,args.allow_underresolved,args.raster)
        elif args.cmd=='run':
            if args.threads<1 or (args.limit_traces is not None and args.limit_traces<1): raise ValueError('Threads and trace limit must be positive')
            run_jobs(args.folder,args.limit_traces,args.threads,args.first_job)
        elif args.cmd=='collect': collect(args.folder)
        elif args.cmd=='grid-sampling': grid_sampling(args.cube,args.out,args.trace_spacings,args.line_spacings,args.gate_start_ns,args.gate_end_ns)
        else: sampling(args.scan,args.out,args.spacings,args.gate_start_ns,args.gate_end_ns)
    except (ValueError,FileNotFoundError) as exc: parser.error(str(exc))
if __name__=='__main__': main()
