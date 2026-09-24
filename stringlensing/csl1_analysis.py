"""Сравнение CSL-1: сдвиговая разность Agol et al. и наклонная струнная модель."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import map_coordinates, maximum_filter
from scipy.optimize import least_squares

import stringlensing as sl


def brightest_pair(image):
    peaks=np.argwhere((image==maximum_filter(image,size=11)) & (image>np.percentile(image,99.4)))
    ordered=sorted(peaks,key=lambda p:image[tuple(p)],reverse=True)
    chosen=[]
    for y,x in ordered:
        if all((x-a)**2+(y-b)**2>15**2 for a,b in chosen): chosen.append((float(x),float(y)))
        if len(chosen)==2:break
    if len(chosen)<2:raise ValueError('Не найдена пара ярких галактик')
    return sorted(chosen)


def pa_moment(image,center):
    x,y=center; yy,xx=np.indices(image.shape)
    rr=np.hypot(xx-x,yy-y)
    mask=(rr>4)&(rr<15)
    sky=np.median(image[rr>20])
    weights=np.clip(image-sky,0,None)*mask
    dx=xx-x;dy=yy-y
    cov=np.array([[np.sum(weights*dx*dx),np.sum(weights*dx*dy)],
                  [np.sum(weights*dx*dy),np.sum(weights*dy*dy)]])/max(np.sum(weights),1)
    values,vectors=np.linalg.eigh(cov);vec=vectors[:,np.argmax(values)]
    return float(np.degrees(np.arctan2(vec[1],vec[0]))%180)


def translation_test(image,noise,first,second,half=15):
    x1,y1=first;x2,y2=second
    gy,gx=np.mgrid[-half:half+1,-half:half+1]
    patch1=map_coordinates(image,[y1+gy,x1+gx],order=1,mode='nearest')
    def second_patch(offset):
        return map_coordinates(image,[y2+offset[1]+gy,x2+offset[0]+gx],order=1,mode='nearest')
    def residual(offset):
        patch2=second_patch(offset)
        variance=2*noise**2+np.clip(patch1,0,None)+np.clip(patch2,0,None)
        return ((patch1-patch2)/np.sqrt(variance)).ravel()
    fit=least_squares(residual,[0.,0.],bounds=([-3,-3],[3,3]),max_nfev=100)
    patch2=second_patch(fit.x)
    chisq=float(np.sum(residual(fit.x)**2))
    return {'translation_dx_pix':float(x2+fit.x[0]-x1),'translation_dy_pix':float(y2+fit.x[1]-y1),
            'redchi':chisq/(patch1.size-2),'chisqr':chisq,'pixels':patch1.size,
            'patch_first':patch1,'patch_second':patch2,'difference':patch1-patch2}


def render_string(p,shape,scale,psf):
    ie,re,n,x0,y0,q,pa,theta,gradient,string_pa,offset,bg=p
    img=sl.generateLensedABSersicGalaxy_new(
        shape[0],shape[1],ie,re,n,x0,y0,q,pa,theta,gradient,string_pa,offset,scale)
    return sl.convolveWithGaussianPSF(img,psf,scale)+bg


def render_independent(p,shape,scale,psf):
    ie1,re1,n1,x1,y1,q1,pa1,ie2,re2,n2,x2,y2,q2,pa2,bg=p
    a=sl.generateSersicGalaxy_new(shape[0],shape[1],ie1,re1,n1,x1,y1,q1,pa1,scale)
    b=sl.generateSersicGalaxy_new(shape[0],shape[1],ie2,re2,n2,x2,y2,q2,pa2,scale)
    return sl.convolveWithGaussianPSF(a+b,psf,scale)+bg


def optimize(image,noise,scale,psf,centers,model,starts,seed,gradient_max,initial_override=None):
    (x1,y1),(x2,y2)=centers
    peak=float(np.max(image));background=float(np.median(image))
    separation=np.hypot(x2-x1,y2-y1)*scale
    angle=float(np.arctan2(x2-x1,-(y2-y1))%np.pi)
    if model=='string':
        initial=np.array([peak/20,.3,2,(x1+x2)/2,(y1+y2)/2,.7,0,separation,0,angle,0,background])
        lo=[.001,.04,.5,0,0,.15,0,.2,-gradient_max,0,-image.shape[1]/2,-peak]
        hi=[peak*3,2.5,6,image.shape[1]-1,image.shape[0]-1,1,math.pi,4,gradient_max,math.pi,image.shape[1]/2,peak]
        render=render_string
    else:
        initial=np.array([peak/20,.3,2,x1,y1,.7,0,peak/20,.3,2,x2,y2,.7,0,background])
        lo=[.001,.04,.5,0,0,.15,0]*2+[-peak]
        hi=[peak*3,2.5,6,image.shape[1]-1,image.shape[0]-1,1,math.pi]*2+[peak]
        render=render_independent
    lo=np.asarray(lo);hi=np.asarray(hi)
    variance=np.maximum(image,0)+noise**2
    valid=np.isfinite(image)
    rng=np.random.default_rng(seed)
    best=None
    for attempt in range(1,starts+1):
        p=np.asarray(initial_override,float).copy() if attempt==1 and initial_override is not None else initial.copy()
        if attempt>1:
            if model=='string':
                p[0]*=rng.uniform(.4,2);p[1]*=rng.uniform(.6,1.6);p[2]=rng.uniform(1,4)
                p[3:5]+=rng.normal(0,3,2);p[6]=rng.uniform(0,math.pi)
                p[7]*=rng.uniform(.75,1.25);p[8]=rng.uniform(-min(300,gradient_max),min(300,gradient_max))
                p[9]=float((angle+rng.normal(0,.3))%math.pi);p[10]=rng.normal(0,5)
            else:
                p[[0,7]]*=rng.uniform(.5,1.8,2);p[[1,8]]*=rng.uniform(.6,1.5,2)
                p[[3,4,10,11]]+=rng.normal(0,2,4)
                p[[6,13]]=rng.uniform(0,math.pi,2)
        p=np.clip(p,lo+1e-7,hi-1e-7)
        def residual(q): return ((render(q,image.shape,scale,psf)-image)/np.sqrt(variance))[valid].ravel()
        try:
            fit=least_squares(residual,p,bounds=(lo,hi),max_nfev=220,ftol=2e-5,xtol=2e-5)
            chisqr=float(np.sum(fit.fun**2))
            if best is None or chisqr<best['chisqr']:
                best={'params':fit.x.tolist(),'chisqr':chisqr,'attempt':attempt,'optimizer_success':bool(fit.success),
                      'nfev':int(fit.nfev),'redchi':chisqr/(valid.sum()-len(p)),
                      'bic':chisqr+len(p)*np.log(valid.sum())}
            print(f'[{model}] start={attempt}/{starts} redchi={chisqr/(valid.sum()-len(p)):.3f} best={best["redchi"]:.3f}',flush=True)
        except Exception as exc:
            print(f'[{model}] start={attempt}: {exc}',flush=True)
    if best is None:raise RuntimeError(f'{model}: все запуски неудачны')
    return best,render(np.asarray(best['params']),image.shape,scale,psf)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-dir',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--starts',type=int,default=8)
    parser.add_argument('--seed',type=int,default=42)
    parser.add_argument('--fit-size',type=int,default=96)
    parser.add_argument('--gradient-max',type=float,default=2000,help='|dThetaE/dXi| в arcsec/rad')
    parser.add_argument('--warm-start',type=Path,help='summary.json предыдущего фита для первого старта струны')
    args=parser.parse_args()
    rows=list(csv.DictReader((args.input_dir/'synthetic_manifest.csv').open()))
    row=next(r for r in rows if json.loads(Path(r['metadata_path']).read_text())['object']=='CSL-1')
    image=np.load(row['image_path']).astype(float)
    pair_full=brightest_pair(image)
    midx=np.mean([p[0] for p in pair_full]); midy=np.mean([p[1] for p in pair_full])
    x0=int(np.clip(round(midx-args.fit_size/2),0,image.shape[1]-args.fit_size))
    y0=int(np.clip(round(midy-args.fit_size/2),0,image.shape[0]-args.fit_size))
    image=image[y0:y0+args.fit_size,x0:x0+args.fit_size]
    noise=float(row['env_noise_sigma']);scale=float(row['env_pixel_scale']);psf=float(row['env_psf_sigma'])
    centers=brightest_pair(image)
    angles=[pa_moment(image,c) for c in centers]
    direct=translation_test(image,noise,*centers)
    warm=json.loads(args.warm_start.read_text())['string']['params'] if args.warm_start else None
    string,string_image=optimize(image,noise,scale,psf,centers,'string',args.starts,args.seed,args.gradient_max,warm)
    independent,independent_image=optimize(image,noise,scale,psf,centers,'independent',args.starts,args.seed+100,args.gradient_max)
    args.output_dir.mkdir(parents=True,exist_ok=True)
    summary={'gradient_bound_arcsec_per_rad':args.gradient_max,
             'warm_start':str(args.warm_start) if args.warm_start else None,
             'source_fits':json.loads(Path(row['metadata_path']).read_text())['source_fits'],
             'centers_xy_pix':centers,'fit_crop_origin_in_prepared_xy':[x0,y0],
             'pixel_scale_arcsec':scale,'psf_fwhm_arcsec':psf,
             'principal_axis_image_deg':angles,'axis_difference_deg':abs((angles[0]-angles[1]+90)%180-90),
             'translation':{k:v for k,v in direct.items() if not isinstance(v,np.ndarray)},
             'string':string,'independent_sersic':independent,
             'delta_bic_string_minus_independent':string['bic']-independent['bic'],
             'caveat':'Оценки chi2 используют диагональную дисперсию; HST drizzle создаёт корреляции и абсолютная значимость не откалибрована.'}
    (args.output_dir/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    fig,axes=plt.subplots(2,3,figsize=(12,8))
    arrays=[image,string_image,image-string_image,independent_image,image-independent_image,direct['difference']]
    titles=['CSL-1 HST F625W','inclined-string Sérsic','string residual',
            'independent Sérsic pair','independent residual','translated difference']
    for ax,arr,title in zip(axes.ravel(),arrays,titles):
        limit=np.percentile(np.abs(arr),99)
        ax.imshow(np.arcsinh(arr/max(noise,1)),origin='lower',cmap='coolwarm' if 'residual' in title or 'difference' in title else 'viridis')
        ax.set_title(title);ax.set_xticks([]);ax.set_yticks([])
    fig.tight_layout();fig.savefig(args.output_dir/'comparison.png',dpi=150);plt.close(fig)
    fig,axes=plt.subplots(1,3,figsize=(10,3.4))
    for ax,key in zip(axes,['patch_first','patch_second','difference']):
        ax.imshow(np.arcsinh(direct[key]/max(noise,1)),origin='lower',cmap='viridis' if key!='difference' else 'coolwarm')
        ax.set_title(key.replace('_',' '));ax.set_xticks([]);ax.set_yticks([])
    fig.tight_layout();fig.savefig(args.output_dir/'translation.png',dpi=150);plt.close(fig)
    print(json.dumps({k:v for k,v in summary.items() if k in {'axis_difference_deg','delta_bic_string_minus_independent'}},indent=2))


if __name__=='__main__':
    main()
