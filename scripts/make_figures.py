from pathlib import Path
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent.parent
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False,'savefig.bbox':'tight'})
d=pd.read_csv(P/'data/v8_28_held_out.csv')
fig,ax=plt.subplots(1,2,figsize=(9,3.4),sharex=True,sharey=True)
for a,c,title in zip(ax,['pred_exchange_RE_Co','pred_nested_selector'],['Fixed RE + Co family','Complete model selection']):
 a.scatter(d.TC_K,d[c],s=30,color='#1b6585',edgecolors='white',linewidth=.4)
 a.plot([30,110],[30,110],'--',color='.4',lw=1);a.set(xlim=(30,110),ylim=(30,110),xlabel='Reported Tc (K)',title=title)
 a.grid(alpha=.18)
 e=d[c]-d.TC_K;a.text(.06,.94,f'MAE = {abs(e).mean():.2f} K\nMax error = {abs(e).max():.2f} K',va='top',transform=a.transAxes,fontsize=9)
ax[0].set_ylabel('Held-out prediction (K)');fig.tight_layout();fig.subplots_adjust(wspace=.20)
for ext in ['png','svg']:fig.savefig(P/f'figures/Figure1_validation.{ext}',dpi=300)
plt.close(fig)
s=pd.read_csv(P/'data/six_frozen_checks.csv');fig,ax=plt.subplots(1,2,figsize=(9,3.4),sharey=True)
for a,el in zip(ax,['Er','Ho']):
 q=s[s.composition.str.contains(el)].copy();q['Gd_pct']=q.composition.str.extract(r'Gd([0-9.]+)').astype(float);q=q.sort_values('Gd_pct')
 a.plot(q.Gd_pct,q.TC_K,'ko-',label='Reported Tc',lw=1)
 a.plot(q.Gd_pct,q.raw_frozen_v5_K,'s--',label='Frozen affine model',color='#b0652b')
 a.plot(q.Gd_pct,q.raw_v8_K,'^--',label='Frozen RE + Co model',color='#1b6585')
 a.set(title=f'Gd–{el} substitution',xlabel='Gd (at.%)',xticks=q.Gd_pct);a.grid(alpha=.18)
ax[0].set_ylabel('Tc (K)');ax[1].legend(fontsize=8,loc='upper left');fig.tight_layout()
for ext in ['png','svg']:fig.savefig(P/f'figures/Figure2_frozen_checks.{ext}',dpi=300)
plt.close(fig)
