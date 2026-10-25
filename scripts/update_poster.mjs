// Replace only the original poster's Figure 1 and Figure 2 image parts.
// Statistical analysis and slide text/geometry are not modified.
import fs from 'node:fs/promises';
import path from 'node:path';
import crypto from 'node:crypto';
import {fileURLToPath} from 'node:url';
import JSZip from 'jszip';

const root=path.resolve(path.dirname(fileURLToPath(import.meta.url)),'..');
const args=process.argv.slice(2);
function option(name,fallback){
  const at=args.indexOf(name);
  if(at<0)return fallback;
  if(!args[at+1])throw new Error('Missing value for '+name);
  return path.resolve(args[at+1]);
}
const source=option('--source',path.join(root,'poster/source/AKI_TAL_original_landscape.pptx'));
const figures=option('--figures',path.join(root,'figures'));
const output=option('--out',path.join(root,'rebuilt/AKI_TAL_poster.pptx'));
if(path.resolve(source)===path.resolve(output))throw new Error('Use a separate output file.');
const digest=b=>crypto.createHash('sha256').update(b).digest('hex');
const before=await fs.readFile(source);
const zip=await JSZip.loadAsync(before);
const replacements=[['ppt/media/image.png','Figure_1_redrawn.png'],['ppt/media/image2.png','Figure_2_redrawn.png']];
for(const [entry,name] of replacements){
  if(!zip.file(entry))throw new Error('Original figure part missing: '+entry);
  const previous=await zip.file(entry).async('nodebuffer');
  const next=await fs.readFile(path.join(figures,name));
  const png=Buffer.from([137,80,78,71,13,10,26,10]);
  if(!previous.subarray(0,8).equals(png)||!next.subarray(0,8).equals(png))throw new Error('Expected PNG figure assets.');
  const oldRatio=previous.readUInt32BE(16)/previous.readUInt32BE(20);
  const newRatio=next.readUInt32BE(16)/next.readUInt32BE(20);
  if(Math.abs(oldRatio-newRatio)>1e-6)throw new Error('Figure aspect ratio must match original: '+name);
  zip.file(entry,next);
}
const bytes=await zip.generateAsync({type:'nodebuffer',compression:'DEFLATE'});
const after=await JSZip.loadAsync(bytes);
const original=await JSZip.loadAsync(before);
const changed=[];
for(const [name,part] of Object.entries(original.files)){
  if(part.dir)continue;
  if(!after.file(name))throw new Error('Package part lost: '+name);
  if(digest(await part.async('nodebuffer'))!==digest(await after.file(name).async('nodebuffer')))changed.push(name);
}
if(JSON.stringify(changed.sort())!==JSON.stringify(replacements.map(r=>r[0]).sort()))throw new Error('Unexpected changed package parts.');
if(digest(await fs.readFile(source))!==digest(before))throw new Error('Source poster changed.');
await fs.mkdir(path.dirname(output),{recursive:true});
await fs.writeFile(output,bytes,{flag:'wx'});
console.log('Poster saved: '+output);
console.log('Updated only Figure 1 and Figure 2. Source poster preserved.');
