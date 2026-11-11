"""General inhumation geometry. SI units; shared definitions with web/cemetery.mjs.
Frames map local x,y,depth points to world x,y,depth. Later graves recut earlier ones.
"""
import copy, json, math
from pathlib import Path
import numpy as np
TEMPLATE=json.loads(Path(__file__).with_name('anatomy.json').read_text())
PLANS=['rectangular','rounded','oval','trapezoid','anthropomorphic']
PROFILES=['vertical','tapered','stepped']
IDENTITY=np.eye(3)

def grave_state(c,g):
 st={k:2**(-g['age_years']/c[k+'_half_life_years']) for k in ['wood','bone','iron','void']}
 if c['environment']=='waterlogged':st['void']=0
 st['collapse']=c['collapse_limit']*(1-st['wood']);return st

def frame(g):
 a=math.radians(g['azimuth_deg']);co,si=math.cos(a),math.sin(a)
 return dict(origin=[g['x_m'],g['y_m'],0],rotation=[co,-si,0,si,co,0,0,0,1])

def rotation_x(a):return np.array([[1,0,0],[0,math.cos(a),-math.sin(a)],[0,math.sin(a),math.cos(a)]])
def rotation_y(a):return np.array([[math.cos(a),0,math.sin(a)],[0,1,0],[-math.sin(a),0,math.cos(a)]])
def world_point(o,p):
 if 'frame' not in o:return list(p)
 f=o['frame'];return (np.array(f['rotation']).reshape(3,3)@p+f['origin']).tolist()

def bounds(o):
 if o['kind']=='loft':
  ps=np.array([[*p,o['z1']] for p in o['lower']]+[[*p,o['z0']] for p in o['upper']]);lo,hi=ps.min(axis=0),ps.max(axis=0)
 elif o['kind']=='box':lo,hi=np.array(o['a']),np.array(o['b'])
 elif o['kind'] in ['sphere','ellipsoid']:lo,hi=np.array(o['p'])-o['r'],np.array(o['p'])+o['r']
 else:lo,hi=np.minimum(o['a'],o['b'])-o['r'],np.maximum(o['a'],o['b'])+o['r']
 if 'frame' not in o:return lo.tolist(),hi.tolist()
 ps=np.array([world_point(o,[hi[i] if k>>i&1 else lo[i] for i in range(3)]) for k in range(8)])
 return ps.min(axis=0).tolist(),ps.max(axis=0).tolist()

def in_polygon(poly,x,y):
 x,y=np.broadcast_arrays(x,y);inside=np.zeros(x.shape,dtype=bool)
 for a,b in zip(poly,poly[-1:]+poly[:-1]):
  if a[1]==b[1]:continue
  inside^=((a[1]>y)!=(b[1]>y))&(x<(b[0]-a[0])*(y-a[1])/(b[1]-a[1])+a[0])
 return inside

def contains(o,X,Y,D):
 X,Y,D=np.broadcast_arrays(X,Y,D)
 if 'frame' in o:
  R=np.array(o['frame']['rotation']).reshape(3,3);p=o['frame']['origin'];dx,dy,dz=X-p[0],Y-p[1],D-p[2]
  X,Y,D=[R[0,i]*dx+R[1,i]*dy+R[2,i]*dz for i in range(3)]
 if o['kind']=='loft':
  t=np.clip((D-o['z0'])/(o['z1']-o['z0']),0,1);inside=np.zeros(X.shape,dtype=bool)
  for i in range(len(o['upper'])):
   j=(i-1)%len(o['upper']);a=[o['upper'][i][k]*(1-t)+o['lower'][i][k]*t for k in range(2)];b=[o['upper'][j][k]*(1-t)+o['lower'][j][k]*t for k in range(2)];den=b[1]-a[1]
   inside^=((a[1]>Y)!=(b[1]>Y))&(X<(b[0]-a[0])*(Y-a[1])/np.where(abs(den)>1e-20,den,1)+a[0])
  return inside&(D>=o['z0'])&(D<o['z1'])
 if o['kind']=='box':return (X>=o['a'][0])&(X<o['b'][0])&(Y>=o['a'][1])&(Y<o['b'][1])&(D>=o['a'][2])&(D<o['b'][2])
 if o['kind']=='ellipsoid':return ((X-o['p'][0])/o['r'][0])**2+((Y-o['p'][1])/o['r'][1])**2+((D-o['p'][2])/o['r'][2])**2<=1
 if o['kind']=='sphere':return (X-o['p'][0])**2+(Y-o['p'][1])**2+(D-o['p'][2])**2<=o['r']**2
 p=np.array(o['a']);v=np.array(o['b'])-p;u=((X-p[0])*v[0]+(Y-p[1])*v[1]+(D-p[2])*v[2])/float(v@v)
 return (u>=0)&(u<=1)&((X-p[0]-u*v[0])**2+(Y-p[1]-u*v[1])**2+(D-p[2]-u*v[2])**2<=o['r']**2)

def footprint(s,inset=0):
 L,W=s['length_m']-2*inset,s['width_m']-2*inset;a,b=L/2,W/2
 if s['plan']=='rectangular':return [[-a,-b],[a,-b],[a,b],[-a,b]]
 if s['plan']=='trapezoid':return [[-a,-b*s['foot_width_ratio']],[a,-b],[a,b],[-a,b*s['foot_width_ratio']]]
 if s['plan']=='anthropomorphic':return [[-a,-b*.48],[a*.38,-b],[a*.68,-b],[a*.68,-b*.48],[a,-b*.48],[a,b*.48],[a*.68,b*.48],[a*.68,b],[a*.38,b],[-a,b*.48]]
 if s['plan']=='oval':return [[a*math.cos(i*math.pi*2/40),b*math.sin(i*math.pi*2/40)] for i in range(40)]
 r=min(a,b);out=[]
 for i in range(34):
  side=1 if i<17 else -1;t=-math.pi/2+(i%17)*math.pi/16+(math.pi if side==-1 else 0)
  out.append([side*(a-r)+r*math.cos(t),r*math.sin(t)])
 return out

def shaft_pieces(s):
 zero=footprint(s);bottom=footprint(s,0 if s['profile']=='vertical' else s['inset_m'])
 if s['profile']=='stepped':
  z=s['top_depth_m']+(s['base_depth_m']-s['top_depth_m'])*s['step_fraction']
  return [dict(z0=s['top_depth_m'],z1=z,upper=zero,lower=zero),dict(z0=z,z1=s['base_depth_m'],upper=bottom,lower=bottom)]
 return [dict(z0=s['top_depth_m'],z1=s['base_depth_m'],upper=zero,lower=bottom)]

def body_objects(g,st,rng):
 if not g['contents']['bones']:return []
 scale=g['body']['stature_m']/1.75;pose=g['body']['pose'];roll=math.pi if pose=='prone' else math.pi/2 if pose.startswith('left') else -math.pi/2 if pose.startswith('right') else 0
 R=rotation_x(roll);out=[];groups={}
 for shape in TEMPLATE:
  name=shape['group']
  if name not in groups:groups[name]=rng.random()<st['bone']
  o=copy.deepcopy(shape);M=IDENTITY;offset=np.zeros(3);leg=name.rsplit('_',1)
  if len(leg)==2 and leg[0] in ['femur','tibia','fibula','patella','foot'] and pose.endswith('flexed'):
   side=int(leg[1]);hip=np.array([.84,side*.083,-.053]);knee=np.array([.48,side*.077,-.047]);A=rotation_y(math.radians(-75))
   if leg[0]=='femur':M=A;offset=hip-A@hip
   else:M=rotation_y(math.radians(65));nk=A@(knee-hip)+hip;offset=nk-M@knee
  for k in ['a','b','p']:
   if k in o:o[k]=((np.array(o[k])-[.90,0,0])*scale).tolist()
  if isinstance(o.get('r'),list):o['r']=[v*scale for v in o['r']]
  elif 'r' in o:o['r']*=scale
  center=np.array([.90,0,0]);shift=(M@center+offset-center)*scale
  o['frame']=dict(rotation=(R@M).ravel().tolist(),origin=(R@shift).tolist());out.append(o)
 if not out:return out
 bs=np.array([bounds(o) for o in out]);lo=bs[:,0,:].min(axis=0);hi=bs[:,1,:].max(axis=0);wall=g['container']['wall_m'] if g['construction']!='earth' else 0
 floor=g['shaft']['base_depth_m']-.035-wall;center=np.array([-(lo[0]+hi[0])/2,-(lo[1]+hi[1])/2,floor-hi[2]])
 for o in out:o['frame']['origin']=(np.array(o['frame']['origin'])+center).tolist()
 return [o for o in out if groups[o['group']]]

def geometry_v3(c,rng_type):
 out=[]
 for i,g in enumerate(c['graves']):
  st=grave_state(c,g);f=frame(g);R=np.array(f['rotation']).reshape(3,3);rng=rng_type(c['seed']+201+g['seed_offset']);own=[]
  def box(a,b,mat,group):own.append(dict(kind='box',a=a,b=b,mat=mat,group=group))
  for p in shaft_pieces(g['shaft']):own.append(dict(kind='loft',**p,mat='fill',group='grave_cut'))
  b=g['container'];L,W,H,t=b['length_m'],b['width_m'],b['height_m'],b['wall_m'];bottom=g['shaft']['base_depth_m']-.035;top=bottom-H
  if g['construction']!='earth':
   wood=g['construction']=='coffin';height=max(H*(1-st['collapse']),g['body']['stature_m']/1.75*(.44 if ('left' in g['body']['pose'] or 'right' in g['body']['pose']) else .24)+2*t) if wood else H;top=bottom-height
   if not wood or st['wood']>.002:
    mat='wood' if wood else 'stone';box([-L/2,-W/2,bottom-t],[L/2,W/2,bottom],mat,'floorboard')
    for side in [-1,1]:
     y=-W/2 if side<0 else W/2-t;box([-L/2,y,top],[L/2,y+t,bottom],mat,'sideboard');x=-L/2 if side<0 else L/2-t;box([x,-W/2,top],[x+t,W/2,bottom],mat,'endboard')
    if g['contents']['lid']:box([-L/2,-W/2,top],[L/2,W/2,top+t],mat,'lid_0')
   if g['contents']['void'] and st['void']>.01 and g['contents']['lid']:
    vh=(bottom-t-top-t)*st['void']
    if vh>.001:box([-L/2+t,-W/2+t,top+t],[L/2-t,W/2-t,top+t+vh],'air','coffin_void')
  own+=body_objects(g,st,rng)
  if g['contents']['iron']:
   z=bottom-t-.025;r=.006*math.sqrt(st['iron']);box([-.025,-.035,z-.012],[.025,.035,z+.012],'rust','iron_object')
   if st['iron']>.005:box([-.018,-.025,z-r],[.018,.025,z+r],'pec','iron_core')
  if g['contents']['pottery']:
   x=0;y=0;z=g['shaft']['top_depth_m']+.10;own+=[dict(kind='sphere',p=[x,y,z],r=.065,mat='pottery',group='vessel'),dict(kind='sphere',p=[x,y,z],r=.057,mat='fill',group='vessel')]
  for o in own:
   if 'frame' in o:
    f2=o['frame'];o['frame']=dict(rotation=(R@np.array(f2['rotation']).reshape(3,3)).ravel().tolist(),origin=(R@f2['origin']+f['origin']).tolist())
   else:o['frame']=copy.deepcopy(f)
   o['grave_id']=g['id']
   if o['mat'] in ['wood','fill']:o['mat']+='_g'+str(i)
   out.append(o)
 return out

def validate_graves(c,rng_type):
 graves=c['graves'];ids=set();d=c['domain']
 if not isinstance(graves,list) or not 1<=len(graves)<=12:raise ValueError('Use 1–12 graves per scene')
 for g in graves:
  def fail(msg):raise ValueError(g.get('label','Grave')+': '+msg)
  s,b=g['shaft'],g['container']
  if not g['id'] or g['id'] in ids:fail('Use unique grave IDs')
  ids.add(g['id'])
  if s['plan'] not in PLANS or s['profile'] not in PROFILES:fail('Unknown shaft shape')
  if g['construction'] not in ['earth','coffin','cist'] or g['body']['pose'] not in ['supine','prone','left_extended','right_extended','left_flexed','right_flexed']:fail('Unknown burial arrangement')
  if g['age_years']<0 or not isinstance(g['seed_offset'],int) or not .45<=g['body']['stature_m']<=2.2:fail('Check age, seed and template stature (0.45–2.2 m)')
  if s['length_m']<=.4 or s['width_m']<=.25 or s['top_depth_m']<0 or s['base_depth_m']<=s['top_depth_m']+.15 or s['base_depth_m']>d['soil_depth_m']-.20:fail('Shaft must fit above the basal margin')
  if s['inset_m']<0 or 2*s['inset_m']>=min(s['length_m']-.2,s['width_m']-.2) or not 0<s['step_fraction']<1 or not .35<=s['foot_width_ratio']<=1:fail('Check taper, bench level or foot width')
  fill=g['fill']
  if not 0<fill['porosity']<1 or not 0<=fill['water_fraction']<=fill['porosity'] or fill['sigma_S_m']<0:fail('Invalid backfill properties')
  for p in shaft_pieces(s):
   lo,hi=bounds(dict(kind='loft',**p,frame=frame(g)))
   if lo[0]<.2 or lo[1]<.2 or hi[0]>d['x_m']-.2 or hi[1]>d['y_m']-.2:fail('Leave 0.2 m around the rotated shaft or enlarge the domain')
  poly=footprint(s,0 if s['profile']=='vertical' else s['inset_m'])
  if g['construction']!='earth':
   if b['length_m']<=0 or b['width_m']<=0 or b['height_m']<=.12 or b['wall_m']<=0 or 2*b['wall_m']>=min(b['length_m'],b['width_m'],b['height_m']):fail('Invalid container dimensions')
   if s['base_depth_m']-.035-b['height_m']<s['top_depth_m']:fail('Container rises above the shaft')
   for x in [-b['length_m']/2,b['length_m']/2]:
    for y in [-b['width_m']/2,b['width_m']/2]:
     if not in_polygon(poly,x*.999999,y*.999999):fail('Container does not fit shaft base')
  bodies=body_objects(g,dict(bone=1),rng_type(1))
  if bodies:
   bs=np.array([bounds(o) for o in bodies]);lo=bs[:,0,:].min(axis=0);hi=bs[:,1,:].max(axis=0)
   if g['construction']!='earth':
    if lo[0]<-b['length_m']/2+b['wall_m'] or hi[0]>b['length_m']/2-b['wall_m'] or lo[1]<-b['width_m']/2+b['wall_m'] or hi[1]>b['width_m']/2-b['wall_m'] or lo[2]<s['base_depth_m']-.035-b['height_m']+b['wall_m']:fail('Body template does not fit container')
   else:
    for x in [lo[0],hi[0]]:
     for y in [lo[1],hi[1]]:
      if not in_polygon(poly,x,y):fail('Body template does not fit shaft base')
   if lo[2]<s['top_depth_m']:fail('Body template rises above shaft')

def voxel_boxes(o,h,ext,surface):
 """Exact cell-centre mask, compressed into non-overlapping axis-aligned boxes.
 Bound memory to one XY plane of one object. Merge equal runs in y, then z.
 """
 lo,hi=bounds(o)
 start=[max(0,int(math.floor(lo[0]/h))),max(0,int(math.floor(lo[1]/h))),max(0,int(math.floor((surface-hi[2])/h)))]
 stop=[min(round(ext[0]/h),int(math.ceil(hi[0]/h))),min(round(ext[1]/h),int(math.ceil(hi[1]/h))),min(round(ext[2]/h),int(math.ceil((surface-lo[2])/h)))]
 if any(stop[k]<=start[k] for k in range(3)):return
 x=(np.arange(start[0],stop[0])+.5)*h;y=(np.arange(start[1],stop[1])+.5)*h;X,Y=np.meshgrid(x,y,indexing='xy');active={}
 for z in range(start[2],stop[2]):
  mask=contains(o,X,Y,surface-(z+.5)*h);rects={};last={}
  for j,row in enumerate(mask):
   edges=np.flatnonzero(np.diff(np.r_[False,row,False].astype(np.int8)));runs=[(int(a)+start[0],int(b)+start[0]) for a,b in zip(edges[::2],edges[1::2])];now={}
   for key in runs:now[key]=last.get(key,(start[1]+j,start[1]+j))[0],start[1]+j+1
   for key,(ya,yb) in last.items():
    if key not in now:rects[(key[0],key[1],ya,yb)]=True
   last=now
  for key,(ya,yb) in last.items():rects[(key[0],key[1],ya,yb)]=True
  for key,za in list(active.items()):
   if key not in rects:
    a,b,ya,yb=key;yield [a*h,ya*h,za*h],[b*h,yb*h,z*h];del active[key]
  for key in rects:active.setdefault(key,z)
 for (a,b,ya,yb),za in active.items():yield [a*h,ya*h,za*h],[b*h,yb*h,stop[2]*h]
