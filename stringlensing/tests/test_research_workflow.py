"""Регрессии координат, физических единиц и вложенного сравнения моделей."""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from photometry_cli import conversion, prepare
from real_fit import components, fit_model, render


class PreprocessingTest(unittest.TestCase):
    def test_hst_rate_uses_exposure_not_ccd_gain(self):
        h=fits.Header();h['BUNIT']='ELECTRONS/S';h['EXPTIME']=120;h['CCDGAIN']=2
        self.assertEqual(conversion(h,None),(120,'electron_rate_times_exposure'))

    def test_wcs_crop_background_gain_and_mask(self):
        rng=np.random.default_rng(7)
        ny=nx=128
        yy,xx=np.indices((ny,nx))
        data=10+rng.normal(0,.7,(ny,nx))+100*np.exp(-((xx-70)**2+(yy-63)**2)/(2*2**2))
        data[67,68]=np.nan
        w=WCS(naxis=2);w.wcs.crpix=[64,64];w.wcs.cdelt=[-.0001,.0001]
        w.wcs.crval=[210.,35.];w.wcs.ctype=['RA---TAN','DEC--TAN']
        ra,dec=w.pixel_to_world_values(70.,63.)
        h=w.to_header();h['GAIN']=2.;h['BUNIT']='ADU'
        with tempfile.TemporaryDirectory() as folder:
            file=Path(folder)/'image.fits';out=Path(folder)/'prepared'
            fits.PrimaryHDU(data.astype('float32'),header=h).writeto(file)
            row,meta=prepare(file,out,'test',float(ra),float(dec),32,psf_fwhm=1.0)
            image=np.load(row['image_path']);mask=np.load(row['valid_mask_path'])
            self.assertEqual(image.shape,(32,32))
            self.assertAlmostEqual(meta['pixel_x'],70.,places=4)
            self.assertAlmostEqual(meta['pixel_y'],63.,places=4)
            self.assertAlmostEqual(meta['conversion_to_electrons'],2.)
            self.assertAlmostEqual(meta['pixel_scale_arcsec'],.36,places=3)
            self.assertGreater(float(image[16,16]),180.)
            self.assertEqual(int(np.count_nonzero(~mask)),1)
            self.assertTrue(np.isfinite(image).all())


    def test_manual_crop_without_wcs_requires_measured_scale(self):
        rng=np.random.default_rng(17)
        yy,xx=np.indices((96,96))
        data=15+rng.normal(0,1,(96,96))+80*np.exp(-((xx-47)**2+(yy-49)**2)/8)
        h=fits.Header();h['BUNIT']='ADU';h['GAIN']=1.7
        with tempfile.TemporaryDirectory() as folder:
            file=Path(folder)/'hct_no_wcs.fits';out=Path(folder)/'prepared'
            fits.PrimaryHDU(data.astype('float32'),header=h).writeto(file)
            with self.assertRaisesRegex(ValueError,'pixel-scale'):
                prepare(file,out,'manual',None,None,32,x=47,y=49,psf_fwhm=1.0)
            row,meta=prepare(file,out,'manual',None,None,32,x=47,y=49,
                             pixel_scale=.3,psf_fwhm=1.0)
            self.assertEqual(np.load(row['image_path']).shape,(32,32))
            self.assertIsNone(meta['wcs_header'])
            self.assertAlmostEqual(meta['pixel_scale_arcsec'],.3)
            self.assertAlmostEqual(meta['conversion_to_electrons'],1.7)


class FitTest(unittest.TestCase):
    def test_free_pair_can_start_from_equal_pair(self):
        rng=np.random.default_rng(3);yy,xx=np.indices((40,40))
        data=render([100,20,20,12,.5,0],'equal_pair',xx,yy,1.7)
        mask=np.ones_like(data,dtype=bool)
        equal=fit_model(data,mask,1.,1.7,'equal_pair',2,rng)
        stars,bg=components(equal['params'],'equal_pair')
        (a,x,y),(b,u,v)=stars
        free=fit_model(data,mask,1.,1.7,'independent_pair',2,rng,[a,b,x,y,u,v,bg])
        self.assertLessEqual(free['chisqr'],equal['chisqr']+1e-6)
        self.assertLess(equal['redchi'],1e-7)

    def test_local_unequal_pair_ignores_remote_bright_source(self):
        rng=np.random.default_rng(11);yy,xx=np.indices((80,80))
        data=render([100,18,38,38,45,44,0],'independent_pair',xx,yy,1.5)
        data+=300*np.exp(-((xx-5)**2+(yy-5)**2)/(2*1.5**2))
        pair=fit_model(data,np.ones_like(data,bool),1.,1.5,'independent_pair',14,rng,
                       center=(39.5,39.5),radius=13,max_sep=14)
        stars,_=components(pair['params'],'independent_pair')
        positions=sorted((x,y) for _,x,y in stars)
        self.assertLess(np.hypot(positions[0][0]-38,positions[0][1]-38),.5)
        self.assertLess(np.hypot(positions[1][0]-45,positions[1][1]-44),.5)
        self.assertEqual(len(pair['trial_chisqr']),14)

    def test_galaxy_component_improves_extended_foreground(self):
        rng=np.random.default_rng(12);yy,xx=np.indices((80,80))
        truth=[120,100,27,49,35,16,0,25,33,40,5]
        data=render(truth,'pair_galaxy',xx,yy,1.1)
        mask=np.ones_like(data,bool)
        pair=fit_model(data,mask,1.,1.1,'independent_pair',5,rng,
                       center=(39.5,39.5),radius=30,max_sep=50)
        galaxy=fit_model(data,mask,1.,1.1,'pair_galaxy',6,rng,
                         center=(39.5,39.5),radius=30,max_sep=50)
        self.assertLess(galaxy['bic'],pair['bic'])


if __name__=='__main__':unittest.main()
