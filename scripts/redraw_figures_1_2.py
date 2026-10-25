#!/usr/bin/env python3
"""Redraw the original poster Figures 1 and 2, retaining all original panels.

All model estimates and UMAP coordinates are read from frozen project tables.
No scientific analysis is rerun. Run with --project PROJECT --out OUTPUT.
"""
import argparse
from pathlib import Path
import json
import numpy as np
import pandas as pd
from PIL import Image
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import PathPatch, Circle, FancyBboxPatch
from matplotlib.path import Path as MPath
from matplotlib.colors import LinearSegmentedColormap

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--project',type=Path,required=True)
p.add_argument('--out',type=Path,required=True)
p.add_argument('--biorender-assets',type=Path,help='Use approved BioRender schematics in panels 1a and 2a')
a=p.parse_args();PROJECT=a.project;OUT=a.out;OUT.mkdir(parents=True,exist_ok=True)
BIO=a.biorender_assets
EXT=PROJECT/'results/extension';RES=PROJECT/'results'
BLUE='#3276A5';RED='#BC4B51';TEAL='#31877C';GOLD='#BD8C39';ORANGE='#C77938'
INK='#20252A';GREY='#727A81';LIGHT='#D7DDE1'
COHORTS=['GSE183276','GSE210622']
GENES=['SLC12A1','UMOD','CLDN16','KCNJ10','PTH1R','CLDN10']
CLASSES=['PT','TAL/LOH','DCT','CNT','CD','POD','EC','Stroma','IMM','Other']
PALETTE=dict(zip(CLASSES,['#629B91',ORANGE,'#8575A7','#C283A2','#99A65B','#C5A344','#568AA9','#9D8972','#727A81','#D2D6D9']))
CMAP=LinearSegmentedColormap.from_list('publication',[BLUE,'#FAFAFA',RED])
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],
 'font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,
 'axes.linewidth':.6,'axes.spines.top':False,'axes.spines.right':False,
 'xtick.major.size':3,'ytick.major.size':3,'xtick.major.width':.6,'ytick.major.width':.6,
 'text.color':INK,'axes.labelcolor':INK,'xtick.color':INK,'ytick.color':INK,
 'pdf.fonttype':42,'svg.fonttype':'none','figure.facecolor':'white','savefig.facecolor':'white'})
def ext(n):return pd.read_csv(EXT/(n+'.csv'))
def result(n):return pd.read_csv(RES/(n+'.csv'))
def cond(dx):return RED if dx=='AKI' else GOLD if dx=='CKD' else BLUE
def tag(fig,letter,x,y):fig.text(x,y,letter,fontsize=12,weight='bold',va='top')
def title(fig,s):fig.text(.055,.98,s,fontsize=13,weight='bold',va='top')
def zero(ax):ax.axvline(0,c='#AEB6BC',lw=.7,ls=(0,(3,3)),zorder=0)
def error(ax,r,y,c,marker='o'):
 ax.errorbar(r.estimate,y,xerr=[[r.estimate-r.ci_low],[r.ci_high-r.estimate]],fmt=marker,
             ms=5,c=c,lw=1.05,capsize=2.7,capthick=.8,mec='white',mew=.5,zorder=3)
def save(fig,name):
 for suffix in ['pdf','svg','png']:fig.savefig(OUT/f'{name}.{suffix}',dpi=400)
 plt.close(fig)
def marker_data(ds,genes):
 m=ext(f'atlas_{ds}_markers');m['z']=m.groupby('gene').mean_expr.transform(lambda v:(v-v.mean())/v.std(ddof=1) if v.std(ddof=1)>0 else v*0)
 return m[m.gene.isin(genes)]
def markers(ax,ds,genes):
 m=marker_data(ds,genes)
 s=ax.scatter(m.broad_class.map(dict(zip(CLASSES,range(10)))),m.gene.map(dict(zip(genes,range(len(genes))))),
       s=m.pct_expr*58+.4,c=m.z,cmap=CMAP,vmin=-2.5,vmax=2.5,edgecolors='#ACB4BA',lw=.18)
 ax.set(xticks=range(10),yticks=range(len(genes)),xlim=(-.6,9.6),ylim=(len(genes)-.5,-.7))
 ax.set_xticklabels(CLASSES,rotation=45,ha='right',fontsize=6.8)
 ax.set_yticklabels(genes,fontstyle='italic',fontsize=6.8);ax.tick_params(length=0)
 ax.grid(axis='y',c='#E8ECEF',lw=.4);ax.set_axisbelow(True)
 for sp in ax.spines.values():sp.set_visible(False)
 ax.set_title(ds,loc='left',pad=8)
 return s
def umap(ax,ds):
 u=ext(f'atlas_{ds}_umap')
 for c in CLASSES:
  d=u[u.broad_class==c];ax.scatter(d.umap1,d.umap2,s=.3,c=PALETTE[c],lw=0,alpha=.7,rasterized=True)
 ax.set_aspect('equal');ax.set_xticks([]);ax.set_yticks([])
 for sp in ax.spines.values():sp.set_visible(False)
 ax.annotate('',(.17,.065),(.04,.065),xycoords='axes fraction',arrowprops={'arrowstyle':'->','lw':.7,'color':INK})
 ax.annotate('',(.04,.22),(.04,.065),xycoords='axes fraction',arrowprops={'arrowstyle':'->','lw':.7,'color':INK})
 ax.text(.08,.005,'UMAP 1',transform=ax.transAxes,fontsize=6)
 ax.text(-.015,.13,'UMAP 2',transform=ax.transAxes,fontsize=6,rotation=90)
 ax.set_title(ds,loc='left',weight='bold',pad=19)
 ax.text(0,1.018,('scRNA-seq' if ds==COHORTS[0] else 'snRNA-seq')+f'   {len(u):,} '+('cells' if ds==COHORTS[0] else 'nuclei')+f'   {u.donor.nunique()} donors',transform=ax.transAxes,fontsize=7.5)

def study_diagram(ax):
 ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
 for x,ds,assay,n,sub,primary in [(.02,COHORTS[0],'scRNA-seq','45 donors / 109,741 cells','12 AKI / 18 reference / 15 CKD','10 AKI / 17 reference'),
   (.54,COHORTS[1],'snRNA-seq','12 donors / 106,971 nuclei','8 AKI / 4 non-AKI','8 AKI / 4 non-AKI')]:
  ax.plot([x,x+.43],[.98,.98],c=TEAL,lw=1.6)
  ax.text(x,.86,('DISCOVERY' if x<.1 else 'INDEPENDENT VALIDATION'),fontsize=8,color=TEAL,weight='bold')
  ax.text(x,.68,f'{ds}   {assay}',fontsize=10,weight='bold')
  ax.text(x,.51,f'Atlas: {n}\n{sub}',fontsize=8,linespacing=1.5)
  ax.text(x,.24,f'Primary TAL: {primary}',fontsize=8,weight='bold')
 ax.text(.02,.04,'Donor pseudobulk: six-gene score, effect and CI, exploratory TAL differential expression',fontsize=8,color=GREY)
 ax.annotate('',(.515,.61),(.47,.61),arrowprops={'arrowstyle':'->','color':GREY,'lw':.9})

def figure1():
 fig=plt.figure(figsize=(11.7,12.0));title(fig,'Figure 1   Study design and human kidney cell-atlas audit')
 tag(fig,'a',.025,.94)
 if BIO:
  ax=fig.add_axes([.18,.625,.64,.335]);ax.imshow(np.asarray(Image.open(BIO/'figure1_biorender.jpg')));ax.axis('off')
 else:
  ax=fig.add_axes([.065,.805,.89,.139]);study_diagram(ax)
  # Tissue support is set apart from the two cell-resolved cohorts.
  ax=fig.add_axes([.065,.755,.89,.041]);ax.axis('off');ax.set(xlim=(0,1),ylim=(0,1))
  ax.plot([0,1],[.99,.99],c=LIGHT,lw=.7)
  ax.text(.0,.66,'CROSS-SECTIONAL TISSUE SUPPORT',color=TEAL,weight='bold',fontsize=7)
  ax.text(.0,.13,'GSE139061: bulk RNA-seq   /   GSE30718: microarray',fontsize=7.5)
  ax.text(.56,.66,'PAIRED TRANSPLANT REPERFUSION SUPPORT',color=TEAL,weight='bold',fontsize=7)
  ax.text(.56,.13,'GSE126805: 41 paired allografts (42 grafts; 163 biopsies overall)',fontsize=7.5)
 for j,ds in enumerate(COHORTS):
  x=.065+j*.47;tag(fig,chr(98+j),x-.04,.604 if BIO else .724);ax=fig.add_axes([x,.396 if BIO else .468,.405,.182 if BIO else .244]);umap(ax,ds)
 fig.legend(handles=[Line2D([],[],marker='o',ls='',ms=4,color=PALETTE[c],label=c) for c in CLASSES],
   loc='center',bbox_to_anchor=(.52,.37 if BIO else .445),ncol=10,frameon=False,fontsize=7,columnspacing=1.1,handletextpad=.4)
 genes=['LRP2','SLC34A1','SLC5A2']+GENES+['SLC12A3','PVALB','CALB1','SCNN1G','ATP6V1B1','SLC4A1','NPHS1','NPHS2','PECAM1','EMCN','DCN','COL1A1','PTPRC']
 tag(fig,'d',.025,.348 if BIO else .415)
 for j,ds in enumerate(COHORTS):
  ax=fig.add_axes([.09+j*.45,.185 if BIO else .224,.365,.145 if BIO else .171]);s=markers(ax,ds,genes)
  if BIO:ax.tick_params(axis='y',labelsize=6)
 cb=fig.colorbar(s,cax=fig.add_axes([.935,.20 if BIO else .249,.009,.113]));cb.set_label('Within-gene z',fontsize=7);cb.ax.tick_params(labelsize=6,length=2)
 fig.text(.09,.14 if BIO else .174,'Dot area: cells detecting gene',fontsize=6.5 if BIO else 7,color=GREY)
 la=fig.add_axes([.29 if BIO else .31,.134 if BIO else .168,.30,.02]);la.axis('off');la.set(xlim=(0,1),ylim=(0,1))
 for x,v in zip([.03,.36,.69],[.25,.5,.75]):la.scatter(x,.5,s=v*58,c=GREY);la.text(x+.045,.5,f'{v:.0%}',va='center',fontsize=7)
 counts=pd.concat([ext(f'atlas_{ds}_counts') for ds in COHORTS]);counts['frac']=counts.n_cells/counts.groupby(['cohort','donor']).n_cells.transform('sum')
 counts['order']=counts.diagnosis.map({'Ref':0,'non_AKI':0,'AKI':1,'CKD':2})
 donor=counts[['cohort','donor','diagnosis','order']].drop_duplicates().sort_values(['cohort','order','donor'])
 key={(r.cohort,r.donor):i for i,r in enumerate(donor.itertuples())};sep=(donor.cohort==COHORTS[0]).sum()
 tag(fig,'e',.025,.124 if BIO else .154);ax=fig.add_axes([.09,.055 if BIO else .062,.48,.058 if BIO else .083]);ax.set_title('Recovered donor composition',loc='left',pad=8 if BIO else 11)
 xx=[key[(r.cohort,r.donor)] for r in counts.itertuples()]
 ax.scatter(xx,counts.broad_class.map(dict(zip(CLASSES,range(10)))),s=counts.frac*60+.2,c=counts.diagnosis.map(cond),lw=0)
 ax.set(yticks=range(10),ylim=(9.6,-.7),xlim=(-.8,len(donor)-.2),xticks=range(len(donor)))
 ax.set_yticklabels(CLASSES,fontsize=5.6 if BIO else 6.5);ax.set_xticklabels(donor.donor,rotation=90,fontsize=4.7);ax.tick_params(length=0)
 ax.axvline(sep-.5,c=LIGHT,lw=.8)
 for sp in ax.spines.values():sp.set_visible(False)
 ax.text(.36,1.025,COHORTS[0],transform=ax.transAxes,ha='center',fontsize=6.5);ax.text(.90,1.025,COHORTS[1],transform=ax.transAxes,ha='center',fontsize=6.5)
 tag(fig,'f',.605,.124 if BIO else .154);ax=fig.add_axes([.655,.055 if BIO else .062,.295,.058 if BIO else .083]);ax.set_title('Broad TAL/LOH recovery',loc='left',pad=8 if BIO else 11)
 q=counts[counts.broad_class=='TAL/LOH'];x=[key[(r.cohort,r.donor)] for r in q.itertuples()]
 ax.scatter(x,np.log1p(q.n_cells),s=12,c=q.diagnosis.map(cond),lw=0)
 for v,ls in [(20,':'),(50,'--')]:ax.axhline(np.log1p(v),c=GREY,lw=.6,ls=ls);ax.text(1.01,np.log1p(v),str(v),transform=ax.get_yaxis_transform(),fontsize=6,va='center',color=GREY)
 ax.set(xlim=(-.8,len(donor)-.2),xticks=range(len(donor)),yticks=np.log1p([20,100,500,2000,10000]),yticklabels=['20','100','500','2,000','10,000'],ylabel='Cells / donor (log1p)')
 ax.set_xticklabels(donor.donor,rotation=90,fontsize=4.1);ax.tick_params(axis='x',length=0);ax.tick_params(axis='y',labelsize=6)
 ax.axvline(sep-.5,c=LIGHT,lw=.8)
 ax.text(.35,1.025,COHORTS[0],transform=ax.transAxes,ha='center',fontsize=6.5);ax.text(.90,1.025,COHORTS[1],transform=ax.transAxes,ha='center',fontsize=6.5)
 fig.legend(handles=[Line2D([],[],marker='o',ls='',ms=3.5,color=c,label=l) for c,l in [(BLUE,'Reference / non-AKI'),(RED,'AKI'),(GOLD,'CKD')]],loc='center',bbox_to_anchor=(.78,.14) if BIO else (.39,.024),ncol=3,frameon=False,fontsize=6.1 if BIO else 6.5)
 fig.text(.065,.009,'Atlas and broad TAL/LOH counts are descriptive. Primary inference uses strict author-annotated TAL; recovered counts do not estimate biological abundance.',fontsize=7,color=GREY)
 save(fig,'Figure_1_redrawn')

def nephron(ax):
 ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
 def path(verts,codes,c,w):ax.add_patch(PathPatch(MPath(verts,codes),facecolor='none',edgecolor=c,lw=w,capstyle='round',joinstyle='round'))
 # Schematic nephron. TAL highlights only the thick ascending segment;
 # descending limb and thin loop bend are grey.
 ax.add_patch(Circle((.10,.79),.052,facecolor='#F6F1ED',edgecolor=GREY,lw=.8))
 t=np.linspace(0,4*np.pi,90);ax.plot(.1+.027*np.cos(t),.79+.024*np.sin(t)*np.cos(t/5),c=GREY,lw=.55)
 pts=[(.15,.79),(.18,.94),(.28,.93),(.25,.80),(.21,.65),(.35,.79),(.34,.60),(.34,.47),(.34,.26),(.34,.18),(.34,.06),(.49,.06),(.49,.18),(.49,.28),(.49,.35),(.49,.43)]
 path(pts,[MPath.MOVETO]+[MPath.CURVE4]*15,'#B5BEC4',4.8)
 path([(.49,.43),(.49,.57),(.49,.72),(.49,.81)],[MPath.MOVETO]+[MPath.CURVE4]*3,ORANGE,6.3)
 path([(.49,.81),(.53,.96),(.60,.70),(.62,.85),(.64,.97),(.70,.92),(.71,.82)],[MPath.MOVETO]+[MPath.CURVE4]*6,'#B5BEC4',4.8)
 ax.plot([.75,.75],[.88,.10],c='#B5BEC4',lw=5,solid_capstyle='round')
 ax.plot([.71,.75],[.82,.74],c='#B5BEC4',lw=4.8,solid_capstyle='round')
 ax.text(.015,.945,'Glomerulus',fontsize=7,color=GREY)
 ax.text(.23,.97,'PT',fontsize=8,color=PALETTE['PT'],weight='bold')
 ax.text(.60,.98,'DCT',fontsize=8,color=PALETTE['DCT'],weight='bold')
 ax.text(.78,.79,'CD',fontsize=8,color=PALETTE['CD'],weight='bold')
 ax.text(.45,.90,'TAL',fontsize=9,color=ORANGE,weight='bold')
 ax.text(.01,.005,'Schematic anatomy; highlighted segment = thick ascending limb',fontsize=6.2,color=GREY)

def programme(ax):
 ax.set(xlim=(0,1),ylim=(0,1));ax.axis('off')
 groups=[('TRANSPORT','SLC12A1   KCNJ10'),('PARACELLULAR','CLDN16   CLDN10'),('DIFFERENTIATION','UMOD   PTH1R')]
 for i,(h,g) in enumerate(groups):
  y=.83-i*.30;ax.text(0,y,h,color=GREY,fontsize=6.5,weight='bold');ax.text(0,y-.14,g,color=ORANGE,fontsize=8.6,fontstyle='italic')

def figure2():
 fig=plt.figure(figsize=(11.7,8.6));title(fig,'Figure 2   Frozen TAL programme across independent cross-sectional cohorts')
 tag(fig,'a',.025,.925)
 if BIO:
  ax=fig.add_axes([.055,.635,.435,.33]);ax.imshow(np.asarray(Image.open(BIO/'figure2_biorender.jpg')));ax.axis('off')
 else:
  fig.text(.065,.92,'Prespecified transport / differentiation programme',fontsize=9,va='top')
  ax=fig.add_axes([.065,.675,.29,.225]);nephron(ax)
  ax=fig.add_axes([.30,.695,.19,.175]);programme(ax)
 tag(fig,'b',.505,.925);fig.text(.55,.92,'Six-gene expression context',fontsize=9,va='top')
 ax=fig.add_axes([.565,.708,.365,.164]);classes=CLASSES[:5]
 mk=pd.concat([marker_data(ds,GENES) for ds in COHORTS]);mk=mk[mk.broad_class.isin(classes)]
 xs=[COHORTS.index(r.cohort)*6+classes.index(r.broad_class) for r in mk.itertuples()]
 ys=mk.gene.map(dict(zip(GENES,range(6))))
 s=ax.scatter(xs,ys,s=mk.pct_expr*95+.4,c=mk.z,cmap=CMAP,vmin=-2,vmax=2,edgecolors='#929CA3',lw=.25)
 ax.set(xticks=[j*6+i for j in range(2) for i in range(5)],yticks=range(6),xlim=(-.7,10.7),ylim=(5.6,-.7))
 ax.set_xticklabels(classes*2,rotation=45,ha='right',fontsize=6.5);ax.set_yticklabels(GENES,fontstyle='italic',fontsize=6.6);ax.tick_params(length=0)
 ax.axvline(5,c=LIGHT,lw=.7)
 ax.text(2,-.92,COHORTS[0],ha='center',fontsize=7);ax.text(8,-.92,COHORTS[1],ha='center',fontsize=7)
 for sp in ax.spines.values():sp.set_visible(False)
 cb=fig.colorbar(s,cax=fig.add_axes([.945,.72,.009,.135]));cb.set_label('Within-gene z',fontsize=6.5);cb.ax.tick_params(labelsize=6,length=2)
 fig.text(.56,.655,'Dot area: detecting cells',fontsize=6.5,color=GREY)
 la=fig.add_axes([.74,.65,.19,.018]);la.axis('off');la.set(xlim=(0,1),ylim=(0,1))
 for x,v in zip([.03,.38,.73],[.25,.5,.75]):la.scatter(x,.5,s=v*95,c=GREY);la.text(x+.06,.5,f'{v:.0%}',fontsize=6.2,va='center')
 ev=result('formal_evidence');ev=ev[(ev.tier=='primary_ge50')&ev.group.isin(['Reference','AKI'])]
 fx=result('formal_effects');fx=fx[(fx.tier=='primary_ge50')&(fx.model=='unadjusted')&(fx.term=='AKI')&fx.comparison.isin(['AKI_vs_Reference','AKI_vs_nonAKI'])]
 for j,ds in enumerate(COHORTS):
  x=.075+j*.235;tag(fig,chr(99+j),x-.05,.597);ax=fig.add_axes([x,.337,.17,.229]);d=ev[ev.dataset==ds]
  vals=[d[d.group==g].TAL_identity for g in ['Reference','AKI']]
  boxes=ax.boxplot(vals,positions=[0,1],widths=.43,patch_artist=True,showfliers=False,
       medianprops={'color':INK,'lw':1},boxprops={'lw':.7},whiskerprops={'lw':.6},capprops={'lw':.6})
  rng=np.random.default_rng(413)
  for k,(v,c) in enumerate(zip(vals,[BLUE,RED])):
   boxes['boxes'][k].set(facecolor=c+'22',edgecolor=c)
   ax.scatter(k+rng.uniform(-.15,.15,len(v)),v,s=17,c=c,edgecolors='white',lw=.4,zorder=3)
  ax.set(xticks=[0,1],xticklabels=[('Reference' if j==0 else 'non-AKI')+f'\n(n = {len(vals[0])})',f'AKI\n(n = {len(vals[1])})'],ylim=(-1.55,1.65),xlim=(-.5,1.5))
  ax.set_title(ds,loc='left',weight='bold',fontsize=8.5,pad=12);ax.tick_params(axis='x',length=0)
  if j==0:ax.set_ylabel('Six-gene TAL score (donor pseudobulk)')
  r=fx[fx.dataset==ds].iloc[0];ax.text(.03,.98,f'Model P = {r.p_value:.2g}',transform=ax.transAxes,va='top',fontsize=7)
 fig.text(.31,.282,'Validation: 8 postmortem AKI\nvs 3 nephrectomy + 1 postmortem non-AKI',fontsize=6.5,color=GREY,va='top')
 tag(fig,'e',.505,.597);fig.text(.55,.59,'Cross-sectional effect estimates',fontsize=9,va='top')
 be=result('bulk_effects');be=be[(be.model=='unadjusted')&(be.analysis_level=='primary')]
 rows=pd.concat([fx,be]);ax=fig.add_axes([.625,.337,.17,.225]);tx=fig.add_axes([.805,.337,.155,.225]);tx.axis('off');tx.set(ylim=(-.7,3.9),xlim=(0,1))
 for i,r in enumerate(rows.itertuples()):
  y=3-i;error(ax,r,y,BLUE if i<2 else TEAL,'s' if i<2 else 'o')
  tx.text(.02,y,f'{r.estimate:.2f} [{r.ci_low:.2f}, {r.ci_high:.2f}]\nP = {r.p_value:.2g}',fontsize=7,va='center',linespacing=1.4)
 zero(ax);ax.axhline(1.5,c=LIGHT,lw=.6);ax.set(yticks=[3,2,1,0],yticklabels=rows.dataset.tolist(),ylim=(-.7,3.9),xlim=(-1.85,.18),xticks=[-1.5,-1,-.5,0],xlabel='AKI effect (95% CI)')
 ax.tick_params(axis='y',length=0,labelsize=7)
 fig.text(.55,.288,'Squares: cell-resolved donors   Circles: tissue-level samples\nCohort-specific models; no pooled effect or pooled P value',fontsize=6.5,color=GREY,va='top')
 tag(fig,'f',.025,.225);fig.text(.075,.221,'Per-gene effect direction and leave-one-donor-out stability',fontsize=9,va='top')
 fg=result('formal_genes');fg=fg[(fg.tier=='primary_ge50')&fg.comparison.isin(['AKI_vs_Reference','AKI_vs_nonAKI'])]
 mat=fg.pivot(index='gene',columns='dataset',values='z_effect').reindex(index=GENES,columns=COHORTS)
 ax=fig.add_axes([.105,.068,.70,.12]);im=ax.imshow(mat,cmap=CMAP,vmin=-2,vmax=2,aspect='auto')
 ax.set(yticks=range(6),yticklabels=GENES,xticks=[0,1],xticklabels=COHORTS);ax.tick_params(length=0);ax.set_yticklabels(GENES,fontstyle='italic',fontsize=7)
 for i in range(6):
  for j in range(2):v=mat.iloc[i,j];ax.text(j,i,f'{v:+.2f}',ha='center',va='center',fontsize=7,color='white' if abs(v)>1.2 else INK)
 cb=fig.colorbar(im,cax=fig.add_axes([.82,.068,.009,.12]));cb.set_label('Gene z-effect',fontsize=7);cb.ax.tick_params(labelsize=6,length=2)
 fig.text(.86,.174,'Negative directions\n6/6 discovery\n5/6 validation\n\nLODO estimates\n100% negative\nin both cohorts',fontsize=7,va='top',linespacing=1.3)
 bg=result('bulk_genes');bg=bg[bg.analysis_level=='primary']
 note=[]
 for ds in ['GSE139061','GSE30718']:
  q=bg[bg.dataset==ds];note.append(f'{ds}: {q.direction.eq("negative").sum()}/{len(q)} genes negative')
 fig.text(.105,.029,'   /   '.join(note),fontsize=7,color=GREY)
 fig.text(.065,.009,'Expression-based programme detectability does not establish biological TAL abundance, TAL selectivity or causality.',fontsize=7,color=GREY)
 save(fig,'Figure_2_redrawn')

if __name__=='__main__':
 figure1();print('Figure 1 complete',flush=True)
 figure2();print('Figure 2 complete',flush=True)
 # Separate editable schematic assets, for use at larger sizes if required.
 for name,draw,size in ([] if BIO else [('Figure_1a_study_design',study_diagram,(10.5,2.0)),('Figure_2a_nephron',nephron,(5.5,4.2))]):
  f=plt.figure(figsize=size);ax=f.add_axes([.05,.07,.9,.86]);draw(ax);f.savefig(OUT/(name+'.svg'));plt.close(f)
