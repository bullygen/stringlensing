"""Локальная многозапусковая подгонка реальных точечных пар и галактики J0146."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from scipy.ndimage import gaussian_filter, maximum_filter
from scipy.optimize import least_squares


def components(p, model):
    if model == 'single':
        a, x, y, bg = p
        return [(a, x, y)], bg
    if model == 'equal_pair':
        a, mx, my, sep, angle, bg = p
        dx = sep/2*np.cos(angle); dy = sep/2*np.sin(angle)
        return [(a, mx-dx, my-dy), (a, mx+dx, my+dy)], bg
    a1, a2, x1, y1, x2, y2, bg = p[:7]
    return [(a1, x1, y1), (a2, x2, y2)], bg


def render(p, model, xx, yy, sigma):
    stars, bg = components(p, model)
    out = np.full(xx.shape, bg, dtype=float)
    for amp, x, y in stars:
        out += amp*np.exp(-((xx-x)**2+(yy-y)**2)/(2*sigma*sigma))
    if model == 'pair_galaxy':
        # Круглый профиль Серсика n=1 (экспоненциальная галактика), свёрнутый с PSF.
        amplitude, gx, gy, re = p[7:11]
        intrinsic = amplitude*np.exp(-1.67834699*np.hypot(xx-gx, yy-gy)/re)
        out += gaussian_filter(intrinsic, sigma=sigma, mode='constant', cval=0.0)
    return out


def peak_guesses(data, sigma, center=None, radius=None):
    ny, nx = data.shape
    cx, cy = center if center is not None else ((nx-1)/2, (ny-1)/2)
    radius = radius if radius is not None else min(nx, ny)/3
    yy, xx = np.indices(data.shape)
    smooth = maximum_filter(data, size=max(3, int(round(sigma*2))))
    candidates = np.argwhere((data == smooth) & np.isfinite(data) &
                             ((xx-cx)**2+(yy-cy)**2 <= radius**2))
    candidates = sorted(candidates, key=lambda q: data[tuple(q)], reverse=True)
    chosen = []
    for y, x in candidates:
        if all((x-qx)**2+(y-qy)**2 > max(1.5, sigma)**2 for qx,qy in chosen):
            chosen.append((float(x), float(y)))
        if len(chosen) == 2: break
    if not chosen: chosen = [(cx, cy)]
    if len(chosen) == 1: chosen.append((chosen[0][0]+max(1.5,sigma), chosen[0][1]))
    return chosen


def fit_model(data, mask, noise, psf_sigma, model, starts, rng, initial_params=None,
              center=None, radius=None, max_sep=None):
    ny, nx = data.shape; yy, xx = np.indices(data.shape)
    cx, cy = center if center is not None else ((nx-1)/2, (ny-1)/2)
    radius = radius if radius is not None else min(nx,ny)/3
    max_sep = max_sep if max_sep is not None else radius*1.5
    valid = mask & np.isfinite(data) & ((xx-cx)**2+(yy-cy)**2 <= radius**2)
    obs = data[valid]
    if len(obs) < 30: raise ValueError('Слишком мало валидных пикселей в локальной области')
    weights = np.maximum(np.sqrt(np.maximum(obs,0)+noise**2), 1e-6)
    bg = float(np.nanmedian(obs))
    (x1,y1),(x2,y2) = peak_guesses(data,psf_sigma,(cx,cy),radius*.85)
    peak = max(float(data[int(round(y1)),int(round(x1))]-bg), 1.0)
    midx,midy = (x1+x2)/2,(y1+y2)/2
    sep = min(max(1.0,float(np.hypot(x2-x1,y2-y1))), max_sep*.9)
    angle = float(np.arctan2(y2-y1,x2-x1))
    lo_x,hi_x = max(1,cx-radius),min(nx-2,cx+radius)
    lo_y,hi_y = max(1,cy-radius),min(ny-2,cy+radius)
    if model == 'single':
        initial=np.array([peak,x1,y1,bg]); lower=[0,lo_x,lo_y,-peak*3];upper=[peak*20,hi_x,hi_y,peak*3]
    elif model == 'equal_pair':
        initial=np.array([peak,midx,midy,sep,angle,bg]);
        lower=[0,lo_x,lo_y,.3,-math.pi,-peak*3];upper=[peak*20,hi_x,hi_y,max_sep,math.pi,peak*3]
    else:
        initial=np.array([peak,peak*.5,x1,y1,x2,y2,bg]);
        lower=[0,0,lo_x,lo_y,lo_x,lo_y,-peak*3];upper=[peak*20,peak*20,hi_x,hi_y,hi_x,hi_y,peak*3]
        if model == 'pair_galaxy':
            initial=np.r_[initial,[peak*.12,cx,cy,4.0]]
            lower += [0,lo_x,lo_y,1.0]; upper += [peak*20,hi_x,hi_y,16.0]
    lower=np.array(lower,float);upper=np.array(upper,float)
    if initial_params is not None: initial=np.asarray(initial_params,float)
    if model == 'pair_galaxy':
        # На HST видны два ядра и слабая промежуточная галактика.
        initial=np.array([peak,peak*.8,27,49,35,16,bg,peak*.12,33,40,4.0])
        lower[2:6]=[21,43,29,10];upper[2:6]=[38,56,44,24]
        lower[8:10]=[25,33];upper[8:10]=[43,46]
    best=None; trials=[]
    for attempt in range(1,starts+1):
        p=initial.copy()
        if attempt > 1:
            if model == 'single':
                p[1:3]+=rng.normal(0,1.5,2)
            elif model == 'equal_pair':
                # Структурированное покрытие малых/средних расстояний и ориентаций.
                p[1:3]=[x1,y1] if attempt%3 else [cx,cy]
                p[3]=min(max_sep*.85,max(.6,psf_sigma)*[.7,1.5,3.0,5.0][(attempt-2)%4])
                p[4]=[-math.pi*.75,-math.pi*.5,-math.pi*.25,0,math.pi*.25,math.pi*.5,math.pi*.75][(attempt-2)%7]
            elif model == 'pair_galaxy':
                p[2:6]+=rng.normal(0,1.2,4)
                p[8:10]+=rng.normal(0,1.5,2)
                p[7]*=rng.uniform(.3,2)
            else:
                if attempt in (2,3,4):
                    p[2:6]=[x1,y1,x2,y2]
                    p[1]=peak*{2:.15,3:.4,4:1.0}[attempt]
                else:
                    a=2*math.pi*((attempt-5)%8)/8
                    d=min(max_sep*.8,max(.8,psf_sigma)*[.8,1.5,2.5,4][(attempt-5)//8%4])
                    p[2:6]=[x1,y1,x1+d*np.cos(a),y1+d*np.sin(a)]
                    p[1]=peak*[.1,.3,.7,1.5][(attempt-5)%4]
            p[0]*=rng.uniform(.7,1.3)
        p=np.clip(p,lower+1e-5,upper-1e-5)
        def resid(q):
            penalty=[]
            if model not in ('single','equal_pair'):
                distance=np.hypot(q[4]-q[2],q[5]-q[3])
                penalty=[max(0,distance-max_sep)*max(peak/noise,1)]
            return np.r_[(render(q,model,xx,yy,psf_sigma)[valid]-obs)/weights,penalty]
        try:
            result=least_squares(resid,p,bounds=(lower,upper),max_nfev=160,ftol=3e-5,xtol=3e-5)
            chi=float(np.sum(result.fun**2));trials.append(round(chi,3))
            if best is None or chi<best['chisqr']:
                best={'params':result.x.tolist(),'chisqr':chi,'attempt':attempt,
                      'optimizer_success':bool(result.success),'message':result.message,
                      'nfev':int(result.nfev)}
        except (ValueError,FloatingPointError):
            continue
    if best is None: raise RuntimeError(f'Все запуски {model} завершились ошибкой')
    n=int(valid.sum());k=len(best['params'])
    best.update(redchi=best['chisqr']/max(1,n-k),bic=best['chisqr']+k*np.log(n),
                n_valid=n, trial_chisqr=trials, distinct_minima=len(set(round(x,1) for x in trials)),
                fit_center_xy=[cx,cy],fit_radius_pixels=radius)
    return best


def fit_one(row,starts,rng,assumed_hst,assumed_ps1):
    meta=json.loads(Path(row['metadata_path']).read_text())
    data=np.load(row['image_path']).astype(float)
    mask=np.load(row['valid_mask_path']).astype(bool) if row.get('valid_mask_path') else np.ones_like(data,bool)
    if data.shape!=(int(row['env_ny']),int(row['env_nx'])): raise ValueError('Размер массива расходится с manifest')
    scale=float(row['env_pixel_scale']);noise=float(row['env_noise_sigma'])
    is_hst='HST' in meta['source_fits'] or 'drc.fits' in meta['source_fits'].lower()
    psf=float(row['env_psf_sigma']) if row.get('env_psf_sigma') else (assumed_hst if is_hst else assumed_ps1)
    psf_status='measured' if row.get('env_psf_sigma') else 'assumed'
    sigma=psf/(2*np.sqrt(2*np.log(2))*scale)
    center=((data.shape[1]-1)/2,(data.shape[0]-1)/2)
    radius=30 if is_hst else 13
    max_sep=50 if is_hst else 14
    kwargs=dict(center=center,radius=radius,max_sep=max_sep)
    results={}
    for model in ['single','equal_pair']:
        results[model]=fit_model(data,mask,noise,sigma,model,2 if model=='single' else starts,rng,**kwargs)
    equal_stars,equal_bg=components(results['equal_pair']['params'],'equal_pair')
    (a1,x1,y1),(a2,x2,y2)=equal_stars
    nested_start=[a1,a2,x1,y1,x2,y2,equal_bg]
    results['independent_pair']=fit_model(data,mask,noise,sigma,'independent_pair',starts,rng,initial_params=nested_start,**kwargs)
    if meta['object']=='J0146–1133' or meta['object']=='J0146-1133':
        results['pair_galaxy']=fit_model(data,mask,noise,sigma,'pair_galaxy',max(6,starts),rng,**kwargs)
    eq=results['equal_pair']; free=results['independent_pair']
    pair,_=components(eq['params'],'equal_pair')
    (a1,x1,y1),(a2,x2,y2)=pair
    sep_arcsec=float(np.hypot(x2-x1,y2-y1)*scale)
    PA_img_deg=float((np.degrees(np.arctan2(y2-y1,x2-x1))+90)%180)
    edge=any(min(x,y,data.shape[1]-1-x,data.shape[0]-1-y)<3 for _,x,y in pair)
    amp_snr=float(a1/max(noise,1e-6))
    row_out={'sample_id':row['sample_id'],'object':meta['object'],'source_fits':meta['source_fits'],
             'psf_fwhm_arcsec':psf,'psf_status':psf_status,'noise_sigma_electrons':noise,
             'equal_redchi':eq['redchi'],'free_redchi':free['redchi'],'single_redchi':results['single']['redchi'],
             'delta_bic_equal_minus_free':eq['bic']-free['bic'],
             'delta_bic_single_minus_free':results['single']['bic']-free['bic'],
             'equal_separation_arcsec':sep_arcsec,'string_angle_image_deg':PA_img_deg,
             'equal_amplitude_electrons':a1,'equal_peak_snr':amp_snr,
             'fit_adequate':bool(eq['optimizer_success'] and eq['redchi']<3 and not edge and amp_snr>5 and sep_arcsec>.25),
             'edge_solution':bool(edge),'n_valid':eq['n_valid'],'best_start':eq['attempt'],
             'equal_distinct_minima':eq['distinct_minima'],'free_distinct_minima':free['distinct_minima'],
             'galaxy_redchi':results.get('pair_galaxy',{}).get('redchi',''),
             'galaxy_delta_bic_pair_minus_galaxy':free['bic']-results['pair_galaxy']['bic'] if 'pair_galaxy' in results else ''}
    return row_out,results


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-dir',type=Path,required=True)
    p.add_argument('--output-dir',type=Path,required=True)
    p.add_argument('--starts',type=int,default=12)
    p.add_argument('--rounds',type=int,default=1)
    p.add_argument('--seed',type=int,default=42)
    p.add_argument('--only-object')
    p.add_argument('--assumed-hst-psf',type=float,default=.10)
    p.add_argument('--assumed-ps1-psf',type=float,default=1.2)
    args=p.parse_args()
    manifest=list(csv.DictReader((args.input_dir/'synthetic_manifest.csv').open()))
    args.output_dir.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(args.seed)
    manifest_path=args.output_dir/'fit_manifest.csv'
    previous={r['sample_id']:r for r in csv.DictReader(manifest_path.open())} if manifest_path.exists() else {}
    payload_dir=args.output_dir/'json';payload_dir.mkdir(exist_ok=True)
    for round_number in range(1,args.rounds+1):
        improved=0
        for row in manifest:
            meta=json.loads(Path(row['metadata_path']).read_text())
            if meta['object']=='CSL-1' or (args.only_object and meta['object']!=args.only_object): continue
            sid=row['sample_id']
            try:
                summary,results=fit_one(row,args.starts,rng,args.assumed_hst_psf,args.assumed_ps1_psf)
                old=previous.get(sid)
                if old is None or float(summary['equal_redchi'])<float(old['equal_redchi']):
                    previous[sid]=summary;improved+=1
                    (payload_dir/f'{sid}.json').write_text(json.dumps(results,ensure_ascii=False,indent=2),encoding='utf-8')
                print(f"[round {round_number}] {sid}: equal χ²ν={summary['equal_redchi']:.2f}; free χ²ν={summary['free_redchi']:.2f}; best={results['equal_pair']['attempt']}/{args.starts}",flush=True)
            except Exception as exc: print(f'[failed] {sid}: {exc}',flush=True)
        if previous:
            with manifest_path.open('w',newline='',encoding='utf-8') as f:
                writer=csv.DictWriter(f,fieldnames=list(next(iter(previous.values()))));writer.writeheader();writer.writerows(previous.values())
        print(f'[round {round_number}] improved={improved}, total={len(previous)}',flush=True)


if __name__=='__main__': main()
