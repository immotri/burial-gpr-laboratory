"""Version 2 geometry using one shared procedural anatomy template."""
import copy,json,math
from pathlib import Path
TEMPLATE=json.loads(Path(__file__).with_name('anatomy.json').read_text())

def geometry_v2(c,st,rng):
 b,d,tar=c['coffin'],c['domain'],c['targets'];xc,yc=d['x_m']/2,d['y_m']/2;L,W,H,t=b['length_m'],b['width_m'],b['height_m'],b['board_m'];x0,x1=xc-L/2,xc+L/2;y0,y1=yc-W/2,yc+W/2;bottom=b['top_depth_m']+H;floor=bottom-t
 minimum=.28 if tar['bones'] else (.24 if tar['pottery'] else 2*t);top=bottom-max(H*(1-st['collapse']),minimum);out=[]
 def box(a,b,mat,group):out.append(dict(kind='box',a=list(a),b=list(b),mat=mat,group=group))
 def cyl(a,b,r,mat,group):out.append(dict(kind='cylinder',a=list(a),b=list(b),r=r,mat=mat,group=group))
 xx=lambda u:x0+.12+u*(L-.24)/1.85
 # Stepped sloping grave walls. Each lower lift narrows towards the coffin.
 for i in range(7):
  za=c['soil']['topsoil_thickness_m']+(bottom+.08-c['soil']['topsoil_thickness_m'])*i/7;zb=c['soil']['topsoil_thickness_m']+(bottom+.08-c['soil']['topsoil_thickness_m'])*(i+1)/7
  margin=.12+.08*(1-i/6);box((x0-margin,y0-margin,za),(x1+margin,y1+margin,zb),'fill','grave_cut')
 if tar['void'] and st['void']>.01:
  vh=max(.003,(floor-top-t)*st['void']);box((x0+t,y0+t,top+t),(x1-t,y1-t,min(floor,top+t+vh)),'air','coffin_void')
 if tar['wood'] and st['wood']>.002:
  # Separate fitted floorboards, side boards and sagging transverse lid planks.
  for j in range(5):box((x0,y0+W*j/5,bottom-t),(x1,y0+W*(j+1)/5-.002,bottom),'wood','floorboard_'+str(j))
  for side in [0,1]:
   for j in range(3):
    ya=y0 if side==0 else y1-t;za=top+(bottom-top)*j/3
    box((x0,ya,za),(x1,ya+t,top+(bottom-top)*(j+1)/3-.002),'wood',f'sideboard_{side}_{j}')
  for i in range(2):
   xa=x0 if i==0 else x1-t
   box((xa,y0,top),(xa+t,y1,bottom),'wood','endboard_'+str(i))
  for i in range(14):
   xa=x0+L*i/14;xb=x0+L*(i+1)/14-.002;drop=min((floor-top-t)*.4,.03)*st['collapse']*math.sin(math.pi*(i+.5)/14)**2
   box((xa,y0,top+drop),(xb,y1,top+t+drop),'wood','lid_'+str(i))
 if tar['bones']:
  groups={}
  for shape in TEMPLATE:
   name=shape['group']
   if name not in groups:groups[name]=rng.random()<st['bone']
   if not groups[name]:continue
   o=copy.deepcopy(shape)
   for key in ['p','a','b']:
    if key in o:o[key]=[xx(o[key][0]),yc+o[key][1],floor+o[key][2]]
   # x semi-axes follow the longitudinal stature scaling.
   if o['kind']=='ellipsoid':o['r'][0]*=(L-.24)/1.85
   out.append(o)
 if tar['iron']:
  sy=yc-.247;bz=floor-.042
  box((xx(.57),sy-.028,bz-.018),(xx(1.49),sy+.028,bz+.018),'rust','sword_corrosion')
  if st['iron']>.005:
   thick=.006*math.sqrt(st['iron']);w=.036*math.sqrt(st['iron'])
   for i in range(12):
    u=.57+i*.066;wid=w*min(1,(i+.5)/3)
    box((xx(u),sy-wid/2,bz-thick/2),(xx(u+.066),sy+wid/2,bz+thick/2),'pec','sword_blade_'+str(i))
   box((xx(1.34),sy-.07,bz-thick/2),(xx(1.365),sy+.07,bz+thick/2),'pec','sword_guard')
   box((xx(1.365),sy-.010,bz-.008),(xx(1.49),sy+.010,bz+.008),'pec','sword_grip')
 if tar['pottery']:
  for n,(u,yy,r,height) in enumerate([(.19,-.195,.064,.15),(.12,.18,.060,.12),(.32,.18,.045,.09)]):
   # Axisymmetric stepped vessel profile, with a closed base and open mouth.
   for k in range(16):
    a=k/16;b=(k+1)/16;z0=floor-height*a;z1=floor-height*b
    radius=r*(.55+.45*math.sin(math.pi*min(1,(a+.03)/.95))) if a<.75 else r*(.35+.12*(a-.75)/.25)
    cyl((xx(u),yc+yy,z1),(xx(u),yc+yy,z0),radius,'pottery',f'vessel_{n}_outer_{k}')
    if k>0:cyl((xx(u),yc+yy,z1-.0001),(xx(u),yc+yy,z0),max(.008,radius-.008),'fill',f'vessel_{n}_inner_{k}')
 return out
