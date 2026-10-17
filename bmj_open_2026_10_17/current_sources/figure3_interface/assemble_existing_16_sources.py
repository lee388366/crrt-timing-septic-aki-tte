from pathlib import Path
import csv,json,hashlib,copy,difflib,shutil
import os
root=Path(os.environ["CRRT_WORKSPACE_ROOT"]);p=Path(os.environ["CRRT_FIG3_INTERFACE_OUT"]);i=p/"NEW_inputs";i.mkdir(parents=True,exist_ok=True)
def read(f):return list(csv.DictReader(Path(f).open()))
def unique(rows,**keys):
 found=[r for r in rows if all(r[k]==v for k,v in keys.items())];assert len(found)==1,(keys,len(found));return found[0]
primfile=Path(os.environ['CH01_PRIMARY_CSV']);prim=read(primfile)
sofafile=Path(os.environ['CRRT_FIG3_SOFA_CSV']);sofa=read(sofafile)
anchorfile=Path(os.environ['CRRT_FIG3_ANCHOR_CSV']);anchors=read(anchorfile)
bridgefile=Path(os.environ['CRRT_FIG3_BRIDGE_CSV']);bridge=read(bridgefile)
pathfile=Path(os.environ['CRRT_FIG3_PATH_CSV']);paths=read(pathfile)
rescuefile=Path(os.environ['CRRT_FIG3_RESCUE_CSV']);rescue=read(rescuefile)
oldfore=Path(os.environ['CRRT_FIG3_FOREST_TEMPLATE']);oldpr=Path(os.environ['CRRT_FIG3_PATH_RESCUE_TEMPLATE']);olde=Path(os.environ['CRRT_FIG3_ERA_TEMPLATE']);lineage={}
def patch(row):
 key=row['analysis'];z=copy.copy(row)
 if key in ['Primary SOFA 8-14 equipoise','Full strictK3_v2 primary cohort']:
  f=primfile;rr=unique(prim,population='Equipoise SOFA 8-14'if key.startswith('Primary')else'Full cohort');cols=['hr','ci_lower','ci_upper','p_value']
 elif key in [r['analysis']for r in sofa]:
  f=sofafile;rr=unique(sofa,analysis=key);cols=['weighted_hr','ci_lower','ci_upper','p_value']
 elif key.startswith('KDIGO'):
  f=anchorfile;rr=unique(anchors,analysis='legacy K3-6h anchor (NEW CH-01)'if key.startswith('KDIGO3')else'KDIGO2 anchor (NEW CH-01)',population='SOFA 8-14');cols=['HR','ci_lower','ci_upper','p_value']
 elif key.startswith('L0-only'):
  f=bridgefile;rr=unique(bridge,analysis='L0-only denominator (NEW CH-01)',population='SOFA 8-14');cols=['HR','ci_lower','ci_upper','p_value']
 elif key=='SOFA 8-14: IPCW includes stage3_path':
  f=pathfile;rr=unique(paths,analysis=key);cols=['HR','CI_lo','CI_hi','p']
 else:
  f=rescuefile;rr=unique(rescue,analysis=key);cols=['weighted_hr','ci_lower','ci_upper','p_value']
 assert float(rr['n_patients'])==float(z['n_patients'])
 if 'early_n'in rr:assert float(rr['early_n'])==float(z['early_n'])
 for target,field in zip(['weighted_hr','ci_lower','ci_upper','p_value'],cols):z[target]=rr[field]
 lineage[key]={'source':str(f),'source_row':{k:rr[k]for k in ['analysis','population']if k in rr},'source_columns':cols,'n_exact':True,'early_n_same_original_key':z['early_n']}
 return z
def write(f,rows):
 with f.open('w',newline='')as o:w=csv.DictWriter(o,fieldnames=list(rows[0]),lineterminator='\n');w.writeheader();w.writerows(rows)
for source,name in [(oldfore,'forest.csv'),(oldpr,'path_rescue.csv')]:
 rows=read(source);out=[patch(r)for r in rows];assert len(rows)==len(out);write(i/name,out)
old=read(olde);assert len(old)==1;z=old[0].copy();rr=unique(bridge,analysis='Calendar-era stratified Cox (NEW CH-01)',population='SOFA 8-14')
for target,field in [('HR','HR'),('CI_lo','ci_lower'),('CI_hi','ci_upper'),('p_value','p_value'),('ph_global_p','ph_global_p')]:z[target]=rr[field]
write(i/'era.csv',[z]);lineage['Calendar-era stratified']={'source':str(bridgefile),'original_source':'09_工作稿/out/run15/run15_NEW_anchors_L0_era.R L82–86; original cox1 AST and CH01 columns','primary_rows_exact':True}
(i/'ch01').mkdir(exist_ok=True);shutil.copyfile(primfile,i/'ch01/ch01_primary_dod.csv')
(p/'16_ROW_SOURCE_LINEAGE.json').write_text(json.dumps({'unique_scientific_rows':len(lineage),'rows':lineage,'no_new_scientific_estimates_or_analyses':'Only original era/L0 bridge cox1/mk calls rerun after OLD4rowbyteexact; all other values read already gated corrected branch outputs.','original_labels_groups_order_and_N_early_preserved':True,'source_sha256':{str(f):hashlib.sha256(f.read_bytes()).hexdigest()for f in [primfile,sofafile,anchorfile,bridgefile,pathfile,rescuefile]}},ensure_ascii=False,indent=2))
