"""Procedural surface/volume template, not a scan or osteometric reconstruction.
Metric coordinates: x along body from feet; y lateral; z depth relative to coffin floor.
The same generated template is read by the browser and the Python FDTD exporter.
"""
import json,math
from pathlib import Path
root=Path(__file__).resolve().parents[1]
shapes=[]
def primitive(group,kind,mat='bone',**kw):shapes.append(dict(group=group,kind=kind,mat=mat,**kw))
def ell(g,p,r,mat='bone',**kw):primitive(g,'ellipsoid',mat,p=list(p),r=list(r),**kw)
def cyl(g,a,b,r,mat='bone',**kw):primitive(g,'cylinder',mat,a=list(a),b=list(b),r=r,**kw)
def curve(g,pts,r,mat='bone'):
 for a,b in zip(pts[:-1],pts[1:]):cyl(g,a,b,r,mat)
def shaft(g,a,b,r):
 cyl(g,a,b,r)
 for p in [a,b]:ell(g,p,[r*1.55,r*1.55,r*1.3])
 v=[b[i]-a[i] for i in range(3)];ln=math.sqrt(sum(q*q for q in v));u=min(.08,.025/ln)
 cyl(g,[a[i]+u*v[i] for i in range(3)],[b[i]-u*v[i] for i in range(3)],max(.003,r-.004),'fill')
# Cranial vault, facial bones and a separate open mandible.
ell('cranium',[1.705,0,-.104],[.102,.076,.086])
ell('cranium',[1.707,0,-.104],[.094,.068,.078],'fill')
ell('cranium',[1.635,0,-.105],[.046,.060,.047])
for side in [-1,1]:
 ell('cranium',[1.653,side*.030,-.140],[.026,.023,.024],'fill',visible_cavity=True)
 ell('cranium',[1.633,side*.051,-.125],[.018,.012,.022])
ell('cranium',[1.622,0,-.136],[.018,.014,.020],'fill',visible_cavity=True)
pts=[(1.611-.018*math.cos(a),.053*math.sin(a),-.101-.018*math.cos(a)) for a in [math.pi/2+i*math.pi/12 for i in range(13)]]
curve('mandible',pts,.007)
for side in [-1,1]:
 cyl('mandible',[1.611,side*.053,-.100],[1.651,side*.053,-.111],.007)
for i in range(10):
 a=math.pi*(i+.5)/10
 ell('mandible',[1.601-.014*math.sin(a),.043*math.cos(a),-.113],[.005,.004,.006])
# Twenty-four vertebral bodies with transverse and spinous processes.
for i in range(24):
 x=.89+i*.0275;radius=.017-(i/23)*.004
 ell(f'vertebra_{i}',[x,0,-.043],[.011,radius,.014])
 for side in [-1,1]:cyl(f'vertebra_{i}',[x,0,-.043],[x-.003,side*(radius+.014),-.040],.004)
# Twelve pairs of curved ribs, tapering caudally; anterior end at elevated sternum.
for i in range(12):
 x=1.42-i*.025;width=.102+.033*math.sin(math.pi*(i+1)/14);width*=min(1, .65+(11-i)*.14)
 for side in [-1,1]:
  pts=[]
  for k in range(15):
   a=math.pi*k/14
   pts.append([x-.023*math.sin(a),side*(.018+width*math.sin(a)), -.043-.099*(1-math.cos(a))/2])
  curve(f'rib_{side}_{i}',pts,.0035 if i>7 else .004)
for i in range(4):ell('sternum',[1.40-i*.033,0,-.145],[.020,.012,.007])
for side in [-1,1]:
 curve('clavicle_'+str(side),[[1.46,0,-.126],[1.47,side*.07,-.120],[1.45,side*.15,-.094]],.007)
 ell('scapula_'+str(side),[1.41,side*.10,-.035],[.064,.042,.009])
 # Humerus, independent radius and ulna, then small carpals and finger rays.
 shoulder=[1.44,side*.155,-.061];elbow=[1.19,side*.196,-.043];wrist=[.96,side*.187,-.041]
 shaft('humerus_'+str(side),shoulder,elbow,.012)
 shaft('radius_'+str(side),[elbow[0],elbow[1]+side*.007,elbow[2]],[wrist[0],wrist[1]+side*.006,wrist[2]],.007)
 shaft('ulna_'+str(side),[elbow[0],elbow[1]-side*.008,elbow[2]],[wrist[0],wrist[1]-side*.008,wrist[2]],.007)
 for i in range(5):
  yy=wrist[1]+side*(i-2)*.009
  ell('hand_'+str(side),[.945,yy,-.040],[.009,.007,.006])
  xend=.865+.018*abs(i-2);curve('hand_'+str(side),[[.94,yy,-.040],[.90,yy+side*.004,-.038],[xend,yy+side*.006,-.035]],.0035)
  for j in range(2):cyl('hand_'+str(side),[xend-j*.016,yy+side*.006,-.035],[xend-(j+1)*.016+.003,yy+side*.006,-.034],.003)
 # Innominate shell and foramina, sacral attachment.
 ell('pelvis_'+str(side),[.89,side*.077,-.048],[.083,.055,.031])
 ell('pelvis_'+str(side),[.895,side*.079,-.052],[.063,.041,.024],'fill')
 curve('pelvis_'+str(side),[[.87,side*.116,-.049],[.80,side*.108,-.050],[.79,side*.032,-.052],[.85,side*.026,-.053]],.010)
 # Femur shafts with separate neck and head; tibia/fibula and patellae.
 hip=[.84,side*.083,-.053];knee=[.48,side*.077,-.047];ankle=[.145,side*.076,-.040]
 shaft('femur_'+str(side),hip,knee,.014)
 cyl('femur_'+str(side),hip,[.862,side*.055,-.052],.012);ell('femur_'+str(side),[.862,side*.055,-.052],[.019,.019,.019])
 shaft('tibia_'+str(side),[knee[0],knee[1]-side*.007,knee[2]],ankle,.011)
 shaft('fibula_'+str(side),[knee[0]-.012,knee[1]+side*.020,knee[2]+.003],[ankle[0],ankle[1]+side*.014,ankle[2]],.006)
 ell('patella_'+str(side),[.48,side*.077,-.067],[.020,.017,.008])
 ell('foot_'+str(side),[.119,side*.077,-.037],[.030,.025,.017])
 for i in range(5):
  yy=side*.077+(i-2)*.010
  curve('foot_'+str(side),[[.116,yy,-.036],[.074,yy+(i-2)*.002,-.032],[.039+.007*abs(i-2),yy+(i-2)*.003,-.029]],.0035)
ell('sacrum',[.875,0,-.047],[.048,.031,.015])
# Keep serialisation deterministic and numerical values readable.
def round_values(o):
 if isinstance(o,float):return round(o,9)
 if isinstance(o,list):return [round_values(x) for x in o]
 if isinstance(o,dict):return {k:round_values(v) for k,v in o.items()}
 return o
text=json.dumps(round_values(shapes),separators=(',',':'))
(root/'toolkit/anatomy.json').write_text(text+'\n')
(root/'web/anatomy-data.mjs').write_text('export default '+text+';\n')
print(len(shapes),'anatomical volume primitives')
