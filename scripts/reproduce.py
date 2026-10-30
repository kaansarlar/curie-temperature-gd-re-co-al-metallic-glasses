"""Reconstruct the archived grouped models and frozen arithmetic without refitting to checks."""
from pathlib import Path
import json,re
import numpy as np
import pandas as pd
from scipy.linalg import lstsq
P=Path(__file__).resolve().parent.parent
def balance(groups):
    groups=np.asarray(groups);u,n=np.unique(groups,return_counts=True)
    counts=dict(zip(u,n))
    return np.array([1/(len(u)*counts[g]) for g in groups])

def split(groups):
    g=np.asarray(groups)
    for label in sorted(set(g)):
        yield np.flatnonzero(g!=label),np.flatnonzero(g==label)

def independent_fit(d,model):
    y=d.TC_K.to_numpy();w=balance(d.validation_group)
    if model['family']=='median':
        ix=np.argsort(y);value=y[ix[np.searchsorted(np.cumsum(w[ix]),.5)]]
        return None,float(value)
    features=model['features'];X=d[features].to_numpy(float)
    if model['family']=='linear_dG_Co':
        mu=w@X;s=np.sqrt(w@((X-mu)**2));s=np.where(s>1e-12,s,1.)
        X=np.column_stack([np.ones(len(X)),X])
        penalty=np.diag(np.r_[0.,np.sqrt(model['alpha'])*s])
    else:
        s=np.sqrt(w@(X*X));s=np.where(s>1e-14,s,1.)
        penalty=np.diag(np.sqrt(model['alpha'])*s)
    a=np.vstack([np.sqrt(w)[:,None]*X,penalty])
    b=np.r_[np.sqrt(w)*y,np.zeros(X.shape[1])]
    beta=lstsq(a,b,lapack_driver='gelsy')[0]
    return beta,None

def ipredict(d,model,fit):
    beta,value=fit
    if beta is None:return np.full(len(d),value)
    X=d[model['features']].to_numpy(float)
    if model['family']=='linear_dG_Co':X=np.column_stack([np.ones(len(X)),X])
    return X@beta

def independent_select(d,templates):
    # Reconstruct inner alpha selection, including all candidates and simplicity.
    fitted={};choices=[]
    for order,(name,template) in enumerate(templates.items()):
        alphas=([None] if name=='median' else [0.] if name=='dg_constant' else
                [.01,.1,1.,10.] if name=='affine' else [0.,.001,.01,.1,1.])
        options=[]
        for alpha in alphas:
            m=dict(template);m['alpha']=alpha;pp=np.empty(len(d))
            for tr,te in split(d.validation_group):
                pp[te]=ipredict(d.iloc[te],m,independent_fit(d.iloc[tr],m))
            loss=float(balance(d.validation_group)@abs(pp-d.TC_K))
            options.append((loss,-alpha if alpha is not None else -np.inf,alpha))
        loss,_,alpha=min(options);m=dict(template);m['alpha']=alpha;fitted[name]=m
        choices.append((loss,m['parameters'],0 if name=='affine' else 1,order,name))
    best=min(x[0] for x in choices)
    chosen=min([x for x in choices if x[0]<=best+1.],key=lambda x:(x[1],x[2],x[3]))[-1]
    return fitted,chosen


def main():
 count=0;mx=0.;summary=[]
 for scope,n,filename in [('primary15',15,'v8_15_held_out.csv'),('expanded28',28,'v8_28_held_out.csv')]:
  d=pd.read_csv(P/'data'/filename);assert len(d)==n
  models=json.loads((P/'models'/scope/'full_models.json').read_text())
  folds=json.loads((P/'models'/scope/'outer_models.json').read_text())
  for fold in folds:
   te=d.record_id.isin(fold['test_records']);tr=~te
   assert not set(d.loc[tr,'doi'])&set(d.loc[te,'doi'])
   assert not set(d.loc[tr,'canonical'])&set(d.loc[te,'canonical'])
   inner,chosen=independent_select(d.loc[tr],models)
   assert inner[chosen]['family']==fold['selected']['family']
   for name,m in fold['models'].items():
    assert inner[name].get('alpha')==m.get('alpha')
    p=ipredict(d.loc[te],m,independent_fit(d.loc[tr],m))
    delta=float(abs(p-d.loc[te,'pred_'+name]).max());mx=max(mx,delta);count+=1;assert delta<1e-8
   p=ipredict(d.loc[te],fold['selected'],independent_fit(d.loc[tr],fold['selected']))
   assert max(abs(p-d.loc[te,'pred_nested_selector']))<1e-8
  for name in list(models)+['nested_selector']:
   err=d['pred_'+name]-d.TC_K
   summary.append({'scope':scope,'model':name,'MAE_K':abs(err).mean(),'group_MAE_K':float(balance(d.validation_group)@abs(err)),'RMSE_K':float(np.sqrt((err**2).mean())),'max_K':abs(err).max(),'negative':int((d['pred_'+name]<0).sum())})
 six=pd.read_csv(P/'data/six_frozen_checks.csv');assert len(six)==6
 G={'Gd':15.75,'Tb':10.5,'Dy':85/12,'Ho':4.5,'Er':2.55,'Y':0};delta=0
 for r in six.to_dict('records'):
  f={e:float(v)/100 for e,v in re.findall(r'([A-Z][a-z]?)([0-9.]+)',r['composition'])};D=sum(f.get(e,0)*g for e,g in G.items());fr=sum(f.get(e,0) for e in G);co=f.get('Co',0);assert fr>0
  a=-48.46044798198254+12.637366749109304*D+145.00506427349472*co
  b=D*(9.60424325448772+1.944423420441906*(D/fr-10)/5+0.9524460124016577*(co-.225)/.05)
  delta=max(delta,abs(a-r['raw_frozen_v5_K']),abs(b-r['raw_v8_K']))
 assert delta<1e-10
 (P/'results').mkdir(exist_ok=True)
 pd.DataFrame(summary).to_csv(P/'results/reproduced_metrics.csv',index=False)
 result={'passed':True,'outer_models_reconstructed':count,'max_model_prediction_difference_K':mx,'frozen_arithmetic_max_difference_K':delta,'new_targets_used_to_fit':False,'scope':'DOI/composition-group validation, not independent laboratory validation'}
 (P/'results/reproduction.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
if __name__=='__main__':main()
