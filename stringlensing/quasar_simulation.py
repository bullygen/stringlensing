"""Monte Carlo восстановления пары точечных квазаров для DOT/Keck."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import binomtest

from real_fit import components, fit_model

ROOT=Path(__file__).resolve().parent.parent


def simulate(config,band,seeing,exposure,magnitude,separation,size,rng,return_image=False):
    filt=config['sky_background']['filters'][band]
    sky_rate=float(filt['sky_electrons_per_second_per_pixel'])
    mu=float(filt['sky_mu_mag_arcsec2'])
    scale=float(config['instrument']['pixel_scale_arcsec_per_pixel'])
    pixel_area=float(config['instrument']['pixel_area_arcsec2_per_pixel'])
    base_fwhm=float(config['psf']['cases'][seeing]['fwhm_arcsec'])
    wavelength_m=(float(filt['central_wavelength_angstrom'])*1e-10 if 'central_wavelength_angstrom' in filt
                  else float(filt['sky_central_wavelength_um'])*1e-6)
    diffraction_fwhm=1.028*wavelength_m/float(config['telescope']['aperture_m'])*206264.806247
    fwhm=math.hypot(base_fwhm,diffraction_fwhm)
    sig=fwhm/(2*math.sqrt(2*math.log(2))*scale)
    total_flux=sky_rate/pixel_area*10**(.4*(mu-magnitude))*exposure
    peak=total_flux/(2*math.pi*sig**2)
    yy,xx=np.indices((size,size));center=(size-1)/2
    center_x=center+rng.uniform(-.5,.5);center_y=center+rng.uniform(-.5,.5)
    angle=rng.uniform(0,math.pi);sep_pix=separation/scale
    dx=sep_pix/2*math.cos(angle);dy=sep_pix/2*math.sin(angle)
    signal=peak*(np.exp(-((xx-center_x-dx)**2+(yy-center_y-dy)**2)/(2*sig**2))+
                 np.exp(-((xx-center_x+dx)**2+(yy-center_y+dy)**2)/(2*sig**2)))
    sky=sky_rate*exposure
    read_noise=4.0
    data=rng.poisson(np.maximum(signal+sky,0)).astype(float)-sky+rng.normal(0,read_noise,(size,size))
    noise=math.sqrt(sky+read_noise**2)
    single=fit_model(data,np.ones_like(data,bool),noise,sig,'single',1,rng)
    equal=fit_model(data,np.ones_like(data,bool),noise,sig,'equal_pair',3,rng)
    equal_stars,equal_bg=components(equal['params'],'equal_pair')
    (a1,x1,y1),(a2,x2,y2)=equal_stars
    free=fit_model(data,np.ones_like(data,bool),noise,sig,'independent_pair',2,rng,
                   initial_params=[a1,a2,x1,y1,x2,y2,equal_bg])
    stars,_=components(equal['params'],'equal_pair')
    (_,x1,y1),(_,x2,y2)=stars
    recovered_sep=np.hypot(x2-x1,y2-y1)*scale
    recovered_angle=math.atan2(y2-y1,x2-x1)%math.pi
    angle_error=abs((recovered_angle-angle+math.pi/2)%math.pi-math.pi/2)
    detected=(single['bic']-free['bic']>10)
    adequate=equal['redchi']<3
    recovered=detected and adequate and abs(recovered_sep-separation)<.2*separation and angle_error<math.radians(15)
    result = {'detected':detected,'adequate':adequate,'recovered':recovered,
              'equal_redchi':equal['redchi'],'separation_error_arcsec':recovered_sep-separation,
              'angle_error_deg':math.degrees(angle_error)}
    if return_image:
        result['example_image'] = data
    return result


def wilson(k,n):
    if n==0:return (np.nan,np.nan)
    ci=binomtest(k,n).proportion_ci(confidence_level=.95,method='wilson')
    return ci.low,ci.high


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir',type=Path,default=ROOT/'result_quasar_simulation')
    p.add_argument('--samples',type=int,default=30)
    p.add_argument('--magnitude',type=float,default=20)
    p.add_argument('--exposures',type=float,nargs='+',default=[30,300,1000])
    p.add_argument('--separations',type=float,nargs='+',default=[.5,1.,2.])
    p.add_argument('--seed',type=int,default=42)
    args=p.parse_args()
    rng=np.random.default_rng(args.seed)
    output=[];args.output_dir.mkdir(parents=True,exist_ok=True)
    choices={'DOT':['site_median_seeing','good_site_seeing'],
             'KECK':['combined_seeing_0p7_plus_instrumental_center','instrumental_center']}
    for telescope,cases in choices.items():
        config=json.loads((ROOT/f'params_{telescope}.json').read_text())
        for band in ['V','R','I']:
            for seeing in cases:
                for exposure in args.exposures:
                    for separation in args.separations:
                        samples=[]
                        for i in range(args.samples):
                            try:
                                samples.append(simulate(config,band,seeing,exposure,args.magnitude,separation,64,rng))
                            except Exception as exc:
                                print(f'[failed] {telescope}/{band}/{seeing}/{exposure}/{separation} #{i}: {exc}',flush=True)
                        n=len(samples)
                        if not n:continue
                        k_detect=sum(r['detected'] for r in samples)
                        k_good=sum(r['adequate'] for r in samples)
                        k_recover=sum(r['recovered'] for r in samples)
                        lo,hi=wilson(k_recover,n)
                        output.append({'telescope':telescope,'filter':band,'psf_case':seeing,
                                       'psf_fwhm_arcsec':math.hypot(
                                           float(config['psf']['cases'][seeing]['fwhm_arcsec']),
                                           1.028*(float(config['sky_background']['filters'][band]['central_wavelength_angstrom'])*1e-10
                                                  if telescope=='DOT' else float(config['sky_background']['filters'][band]['sky_central_wavelength_um'])*1e-6)
                                           /float(config['telescope']['aperture_m'])*206264.806247),
                                       'pixel_scale_arcsec':config['instrument']['pixel_scale_arcsec_per_pixel'],
                                       'exposure_s':exposure,'separation_arcsec':separation,'magnitude_each':args.magnitude,
                                       'n':n,'n_detected':k_detect,'n_adequate':k_good,'n_recovered':k_recover,
                                       'recovery_fraction':k_recover/n,'recovery_ci95_low':lo,'recovery_ci95_high':hi,
                                       'median_redchi':float(np.median([r['equal_redchi'] for r in samples])),
                                       'median_abs_separation_error_arcsec':float(np.median([abs(r['separation_error_arcsec']) for r in samples]))})
                        print(f'[point] {telescope}/{band}/{seeing} t={exposure:g}s sep={separation:g}" recovery={k_recover}/{n}',flush=True)
    with (args.output_dir/'simulation_summary.csv').open('w',newline='',encoding='utf-8') as f:
        writer=csv.DictWriter(f,fieldnames=output[0].keys());writer.writeheader();writer.writerows(output)
    (args.output_dir/'run_metadata.json').write_text(json.dumps({'samples':args.samples,'magnitude':args.magnitude,
                'exposures':args.exposures,'separations':args.separations,'seed':args.seed,
                'source_model':'two identical Gaussian PSF point sources',
                'diffraction':'FWHM_eff^2 = FWHM_case^2 + (1.028 lambda/D rad * 206264.806)^2',
                'selection':'delta BIC single-free > 10, equal-pair redchi < 3, separation error < 20%, PA error < 15 deg',
                'limitations':'Gaussian PSF, no field contaminants, read noise 4 e-, same magnitude in all filters; no diffraction wings or intrinsic quasar variability.'},indent=2),encoding='utf-8')
    fig,axs=plt.subplots(2,3,figsize=(12,7),sharey=True)
    for ax,telescope,band in [(axs[i,j],t,['V','R','I'][j]) for i,t in enumerate(['DOT','KECK']) for j in range(3)]:
        for case in choices[telescope]:
            rows=[r for r in output if r['telescope']==telescope and r['filter']==band and r['psf_case']==case and r['exposure_s']==300]
            rows=sorted(rows,key=lambda r:r['separation_arcsec'])
            ax.errorbar([r['separation_arcsec'] for r in rows],[r['recovery_fraction'] for r in rows],
                        yerr=[[r['recovery_fraction']-r['recovery_ci95_low'] for r in rows],
                              [r['recovery_ci95_high']-r['recovery_fraction'] for r in rows]],
                        label=case.replace('_',' '),marker='o',capsize=2)
        ax.set_title(f'{telescope} {band}, 300 s');ax.set_xlabel('separation, arcsec');ax.grid(alpha=.2)
        ax.set_ylim(-.05,1.05)
        if band=='V':ax.set_ylabel('recovery fraction (95% Wilson CI)')
        if telescope=='DOT':ax.legend(fontsize=6)
    fig.tight_layout();fig.savefig(args.output_dir/'recovery_vs_separation.png',dpi=150);plt.close(fig)
    print(f'[summary] {len(output)} parameter points, {sum(r["n"] for r in output)} Monte Carlo images')


if __name__=='__main__':
    main()
