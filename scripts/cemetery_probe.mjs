import fs from 'node:fs';
import {geometry,estimate,spectra,sigmas} from '../web/core.mjs';
import {sectionRaster} from '../web/explorer.mjs';
import {validateGraves} from '../web/cemetery.mjs';
let models=JSON.parse(fs.readFileSync(0,'utf8')),out=[];
for(let c of models){validateGraves(c);let objects=geometry(c),slices=[];for(let plane of ['xy','xz','yz']){let d=c.domain,extent=plane==='xy'?[0,d.x_m,0,d.y_m]:plane==='xz'?[0,d.x_m,0,d.soil_depth_m]:[0,d.y_m,0,d.soil_depth_m],r=sectionRaster(c,objects,null,plane,c.sections,extent,97,73);slices.push({plane,extent,data:Array.from(r.data),names:r.names});}out.push({objects,estimate:estimate(c),spectra:spectra(c,500e6),sigmas:sigmas(c),slices});}
console.log(JSON.stringify(out));
