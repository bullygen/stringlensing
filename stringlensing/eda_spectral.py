"""Просмотр локальных спектральных FITS и сохранение обзорного графика."""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

ROOT=Path(__file__).resolve().parent.parent


def load_spectrum(path: Path):
    with fits.open(path,memmap=False,output_verify='silentfix') as hdul:
        hdul.verify('silentfix')
        for index,hdu in enumerate(hdul):
            if hdu.data is None or not np.issubdtype(np.asarray(hdu.data).dtype,np.number):
                continue
            data=np.asarray(hdu.data,dtype=float)
            while data.ndim>1:
                data=data[0]
            if data.ndim!=1 or data.size<3:
                continue
            hdr=hdu.header
            crval=hdr.get('CRVAL1');cdelt=hdr.get('CDELT1',hdr.get('CD1_1'))
            crpix=float(hdr.get('CRPIX1',1))
            if crval is not None and cdelt is not None:
                axis=float(crval)+(np.arange(data.size)+1-crpix)*float(cdelt)
                unit=str(hdr.get('CUNIT1','согласно FITS'))
            else:
                axis=np.arange(data.size)
                unit='пиксель'
            return axis,data,unit,index
    raise ValueError(f'{path}: спектральный массив не найден')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('files',type=Path,nargs='+',help='Локальные спектральные FITS')
    p.add_argument('--output',type=Path,required=True,help='PNG результат')
    args=p.parse_args()
    fig,axes=plt.subplots(len(args.files),1,figsize=(10,max(3,2.7*len(args.files))),squeeze=False)
    for ax,path in zip(axes[:,0],args.files):
        axis,data,unit,hdu=load_spectrum(path)
        ax.plot(axis,data,lw=.7)
        ax.set_title(f'{path.name}, HDU {hdu}')
        ax.set_xlabel(unit);ax.set_ylabel('Сигнал, единицы FITS')
        ax.grid(alpha=.2)
        print(f'[spectrum] {path}: HDU={hdu}, points={data.size}, finite={np.isfinite(data).sum()}')
    fig.tight_layout();args.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(args.output,dpi=150);plt.close(fig)
    print(f'[saved] {args.output}')


if __name__=='__main__':main()
