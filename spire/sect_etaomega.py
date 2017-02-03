from scipy import ndimage as ndimage
from scipy import interpolate as interpolate
from mpfit import mpfit
from scipy.io.idl import readsav
from math import *
from numpy import *
from astropy import wcs
from matplotlib.cbook import is_numlike
import pyfits,os
import pdb

sqasc2sr=(pi/180./3600.)**2.
jytomks=1e-26
deg2ra=(1./360.)*2.*math.pi
h_planck=6.626e-34
k_boltzmann=1.38e-23
c=3e8

def make_dummy_wcs(ra,dec):
    
    dummy=wcs.WCS()
    dummy.wcs.crpix=array([128.,128.])
    dummy.wcs.crval=array([ra,dec])
    dummy.wcs.cdelt=array([-1,1])*1./3600.
    pa=radians(0.)
    cpa=cos(pa)
    spa=sin(pa)
    dummy.wcs.pc=array([[-cpa,-spa],[-spa,cpa]])
    
    return dummy

def planck_nu(t,nu):

    planck_b=2.*h_planck*nu**3./c**2./(exp(h_planck*nu/k_boltzmann/t)-1.)

    return planck_b

# wavelength: in units of micron!
def herschel_em(wavelength):

    em=0.0336*(wavelength)**(-0.5)+0.273*(wavelength)**(-1.)

    return em

def herschel_temperature(tm1,tm2,freq_ghz):

    wavelength=3e5/freq_ghz
    em=herschel_em(wavelength)
    temperature=(1.-em)*em*tm1+em*tm2

    return temperature

def herschel_itel(tm1,tm2,freq_ghz):

    wavelength=3e5/freq_ghz
    em=herschel_em(wavelength)
    itel=(1.-em)*em*planck_nu(tm1,freq_ghz)+em*planck_nu(tm2,freq_ghz)

    return itel
    
# r50: arcsec
def sim_gauss(dfwhm,ctrxy):

    rsig=dfwhm/(2.*sqrt(2.*log(2.)))
    img_gauss=zeros([257,257])
    imgsize=img_gauss.shape

    for i in range(imgsize[0]):
        for j in range(imgsize[1]):
            img_gauss[i,j]=exp(-(1.*(ctrxy[0]-i)**2.+1.*(ctrxy[1]-j)**2.)/2./rsig**2.)

    return img_gauss

def sim_exp(dfwhm,ctrxy):

    rex=dfwhm/2./log(2.)
    img_exp=zeros([257,257])
    imgsize=img_exp.shape

    for i in range(imgsize[0]):
        for j in range(imgsize[1]):
            img_exp[i,j]=exp(-sqrt(1.*(ctrxy[0]-i)**2.+1.*(ctrxy[1]-j)**2.)/rex)

    return img_exp

def sim_planet(dlight,ctrxy):

    img_planet=zeros([257,257])
    imgsize=img_planet.shape

    for i in range(imgsize[0]):
        for j in range(imgsize[1]):
            if (i-ctrxy[1])**2.+(j-ctrxy[0])**2. <= int(ceil(dlight/2.))**2.:
                img_planet[i,j]=1.
                

    return img_planet

def fit_likelihood(p,fjac=None,x=None,y=None,err=None,rec=None):
    ymod=exp(-1.*((x-p[0])/p[1])**2./2.)
    status=0
    return ([status,(y-ymod)/err])

def fit_gaussian(x,a,b):
    y=exp(-1.*((x-a)/b)**2./2.)
    return y

def save_rsrf_fits():
    
    bol_vres=pyfits.HDUList()
    telrsrf=pyfits.HDUList()
    ntrsrf=rtel[0].header['DSETS___']

    for j in range(2):
        if j == 0:
            ni='SSW'
            deti='SSWD4'
            rtel_copy=rtel[21].copy()
            vres=rtel[21].copy()
        if j == 1:
            ni='SLW'
            deti='SLWC3'
            rtel_copy=rtel[2].copy()
            vres=rtel[2].copy()
        rtel_tmp=rtel_copy.copy()
        vres.update_ext_name(deti)
        rtel_tmp.update_ext_name(deti)
        rtel_tmp.data['rsrf'][:]=double(0)
        rtel_tmp.data['error'][:]=double(0)
        vres.data['rsrf'][:]=double(0)
        vres.data['error'][:]=double(0)
        for i in range(ntrsrf):
            dset=str(i+1).zfill(4)
            ds_last='DS_'+str(rtel[dset].header['DSETS___']-1)
            drange=[rtel[dset].header['DS_0'],rtel[dset].header[ds_last]]
            hdu_temp=rtel[drange[0]:drange[1]]
            hdu=hdu_temp[deti].copy()
            rtel_tmp.data['rsrf'][:]=rtel_tmp.data['rsrf'].copy()+ \
                hdu.data['rsrf']/hdu.data['error']**2.
            rtel_tmp.data['error'][:]=rtel_tmp.data['error']+1./hdu.data['error']**2.

        rtel_tmp.data['rsrf'][:]=rtel_tmp.data['rsrf']/rtel_tmp.data['error']
        rtel_tmp.data['error'][:]=sqrt(1./rtel_tmp.data['error'])
        rtel_tmp.header['EXTNAME']=deti
        rtel_tmp.header['META_0']=deti
        telrsrf.append(rtel_tmp)
        vres.header['EXTNAME']=deti
        vres.header['META_0']=deti
        vres.columns[1].name='VolResp'
        vres.columns[1].unit='V W-1 m2'
        vres.columns[2].name='VolRespError'
        vres.columns[2].unit='V W-1 m2'
        vres.header['TTYPE2']='VolResp'
        vres.header['TUNIT2']='V W-1 m2'
        vres.header['TTYPE3']='VolRespError'
        vres.header['TUNIT3']='V W-1 m2'
        vres.data['rsrf'][:]=rtel_tmp.data['rsrf']*1e-9/(cps[deti].data['pointConv']*1e-26)
        vres.data['error'][:]=sqrt( \
            (rtel_tmp.data['error']*1e-9/(cps[deti].data['pointConv']*1e-26))**2.+
            (rtel_tmp.data['rsrf']*1e-9*cps[deti].data['pointConv']*1e-26/ \
                 (cps[deti].data['pointConv']*1e-26))**2.)
        bol_vres.append(vres)
        rtel_tmp=None
        vres=None

    telrsrf.writeto(plotdir+'telersrf_avg.fits',clobber=True, \
                        output_verify='ignore')
    bol_vres.writeto(plotdir+'bol_vresp.fits',clobber=True, \
                         output_verify='ignore')

    return telrsrf,bol_vres

def find_best_d(ftsfile,d_plt_ini,center=None,caltype=None,obj=None):
# ctrxy: array([x,y]) - pixel location of the object center

    ssw=ftsfile['SSWD4']
    slw=ftsfile['SLWC3']
    sswidx=where(ssw.data['wave'] < fmax_slw)
    sswidx=sswidx[0][1:]
    ssw=ssw.data[sswidx]
    ssw=ssw
    slwidx=where(slw.data['wave'] > fmin_ssw)
    slwidx=slwidx[0][0:-1]
    slw=slw.data[slwidx]

    if not is_numlike(center):
        ctrxy=array([128,128])
    else:
        hdr=ftsfile['SLWC3'].header
        dummywcs=make_dummy_wcs(hdr['RA'],hdr['DEC'])
        xyradec=array([center])
        ctrxy_tmp=dummywcs.wcs_world2pix(xyradec,0)
        ctrxy=copy(ctrxy_tmp[0])

    xs=ssw['wave']
    xl=slw['wave']

    if caltype == 'point':
        ssw_cal=zeros(len(cps['SSWD4'].data['pointConv']))+1.
        slw_cal=zeros(len(cps['SLWC3'].data['pointConv']))+1.
        
    if caltype == 'extended':
        ssw_cal=cps['SSWD4'].data['pointConv'].copy()
        slw_cal=cps['SLWC3'].data['pointConv'].copy()

    ssw_cal=ssw_cal[sswidx]
    slw_cal=slw_cal[slwidx]

    if (obj == 'm83' or obj == 'M83') or \
            (obj == 'm82' or obj == 'M82') or \
            (obj == 'lmc-n159' or obj == 'LMC-N159'):
        scal=cps['SSWD4'].data['pointConv'].copy()
        lcal=cps['SLWC3'].data['pointConv'].copy()
        scal=scal[sswidx]
        lcal=lcal[slwidx]
        ssw['flux'][:]=ssw['flux'].copy()*scal
        ssw['error'][:]=ssw['error'].copy()*scal
        slw['flux'][:]=slw['flux'].copy()*lcal
        slw['error'][:]=slw['error'].copy()*lcal
        sim_func=sim_gauss
    else:
        sim_func=sim_planet


    ssw_wnidx=(wn_ssw[0].data*30. < fmax_slw+5.) & \
        (wn_ssw[0].data*30. > fmin_ssw-5.)
    slw_wnidx=(wn_slw[0].data*30. < fmax_slw+5.) & \
        (wn_slw[0].data*30. > fmin_ssw-5.)
    slw_wnidx[-7]=False
    wn_sfit=wn_ssw[0].data[ssw_wnidx]
    wn_lfit=wn_slw[0].data[slw_wnidx]

    ssw_b=beam_ssw[0].data[ssw_wnidx,:,:].copy()
    slw_b=beam_slw[0].data[slw_wnidx,:,:].copy()
    sum_ssw_b=sum(sum(ssw_b,axis=1),axis=1)
    sum_slw_b=sum(sum(slw_b,axis=1),axis=1)

    refidx=1e9
    chiparam=[]
    chierr=[]
    d_plt_out=None
    d_plt=arange(26)*2.
    dplt=d_plt_ini+d_plt-20.
    dplt=dplt[where(dplt >= 0.)]
    for di in range(len(dplt)):
        d_input=dplt[di]
        planet_mod=sim_func(d_input,ctrxy)*img_mask
        planet_mod=planet_mod.copy()/planet_mod.max()
        planet_area=sum(planet_mod)
        sum_ssw_bp=[]
        sum_slw_bp=[]
        for bi in range(len(sum_ssw_b)):
            sum_ssw_bp.append(sum(ssw_b[bi,:,:]*planet_mod))
        for bi in range(len(sum_slw_b)):
            sum_slw_bp.append(sum(slw_b[bi,:,:]*planet_mod))
        sum_ssw_bp=array(sum_ssw_bp)
        sum_slw_bp=array(sum_slw_bp)
        f_sum_sbp=interpolate.interp1d(wn_sfit*30.,planet_area/sum_ssw_bp, \
                                           bounds_error=False, kind=3, \
                                           fill_value=planet_area/sum_ssw_bp[0])
        f_sum_lbp=interpolate.interp1d(wn_lfit*30.,planet_area/sum_slw_bp, \
                                           bounds_error=False, kind=3, \
                                           fill_value=planet_area/sum_slw_bp[0])
        ssw_corrf=ssw['flux']*ssw_cal*f_sum_sbp(xs)
        slw_corrf=slw['flux']*slw_cal*f_sum_lbp(xl)
        param=sum((slw_corrf-ssw_corrf)**2./100./ \
                      2./((slw['error']*slw_cal*f_sum_lbp(xl))**2.+ \
                              (ssw['error']*ssw_cal*f_sum_sbp(xs))**2.))
        err=sqrt(sum((slw['error']*slw_cal*f_sum_lbp(xl))**2.+ \
                    (ssw['error']*ssw_cal*f_sum_sbp(xs))**2.))
        chiparam.append(param)
        chierr.append(err)
        if param < refidx:
            refidx=param

    chiparam=array(chiparam)
    chierr=array(chierr)
    chiidx=(chiparam == chiparam.min())

    return chiparam,chierr,dplt

def plot_cps():

    fwhm_uranus=1.7
    ctrxy=array([128,128])
    uranus=sim_planet(fwhm_uranus,ctrxy)

    beam_area=[]
    beam_uranus=[]
    
    uranus_area=sum(uranus)
    
    figure(0,figsize=[8.,5.])
    subplot(111)
    title('FTS Telescope RSRF',fontsize='large')

    figure(1,figsize=[8.,8.])
    subplot(111)
    hdu=pyfits.PrimaryHDU()
    hdulist=pyfits.HDUList(hdu)
    for i in range(2):
        if i == 0:
            wn=wn_ssw
            beam=beam_ssw
            fwhm=fwhm_ssw
            iro='blue'
            ni='SSW'
        if i == 1:
            wn=wn_slw
            beam=beam_slw
            fwhm=fwhm_slw
            iro='red'
            ni='SLW'

        detector=side_det[i]
        beam_area=[]
        beam_uranus=[]
        norm=[]
        for j in range(len(wn[0].data)):

            beam_area.append(sum(beam[0].data[j,:,:]))
            beam_uranus.append(sum(beam[0].data[j,:,:]*uranus))

        a_beam=None
        ps_coup=None
        eta_mb=None
        a_beam=array(beam_area)
        ps_coup=array(beam_uranus)
        
        x=wn[0].data*30.
        y=a_beam*sqasc2sr
        figure(0)
        plot(x,y,'--',color=iro,label=ni+' model',linewidth=3.)
        plot(x,1.1*y,'--',color=iro,linewidth=1.)
        plot(x,0.9*y,'--',color=iro,linewidth=1.)
        fill_between(x,1.1*y,0.9*y,color=iro,alpha=0.2)
        for idet in range(len(detector)):
            deti=detector[idet]
            figure(0)
            plot(cps[deti].data['wave'],cps[deti].data['pointConv']*1e-26, \
                     label=deti,linewidth=1.)

#        fill_between(cps[deti].data['wave'], \
#                         (cps[deti].data['pointConv']+ \
#                             cps[deti].data['pointConvError'])*1e-26, \
#                        (cps[deti].data['pointConv']- \
#                              cps[deti].data['pointConvError'])*1e-26, \
#                         color=iro,alpha=0.2)
            fy=interpolate.interp1d(x,y,bounds_error=False,fill_value=y[0])
            gamma=fy(cps[deti].data['wave'])/cps[deti].data['pointConv']/1e-26

            figure(1)
#            if i == 0:
#                plot(cps[deti].data['wave'],gamma,color='grey',label='Dark Sky')
#            else:
            plot(cps[deti].data['wave'],gamma,label=deti)
            if deti == 'SLWC3' or deti == 'SSWD4':
                print deti
                c1=pyfits.Column(name='wave',format='D', \
                                     unit='GHz',array=cps[deti].data['wave'])
                c2=pyfits.Column(name='gamma',format='D', \
                                     unit='--',array=gamma)
                tmptable=pyfits.new_table([c1,c2])
                tmptable.name=ni
                hdulist.append(tmptable)

    hdulist.writeto('../results/Beams/gamma.fits',clobber=True)
    figure(0)
#    pdb.set_trace()
    xlabel('Frequency (GHz)',fontsize='large')
    ylabel('C$_{PS}$ (sr)',fontsize='large')
    xlim([400.,1600.])
    ylim([0.,1.5e-7])
    legend(loc=1,fontsize='small',ncol=3)
    savefig(filename=plotdir+'cps_v11.png',format='PNG')
    os.system('xview -shrink %scps_v11.png &' %plotdir)
    close(0)

    figure(1)

#    bms=readsav('../data/FTS/sect_planets/stuff_needed_SECT_eff_plot.sav')
#    plot(29.98*bms.wns[bms.wsx],bms.yfs,'-.',color='red',linewidth=3.)
#    plot(29.98*bms.wnl[bms.wlx],bms.yfl,'-.',color='red',linewidth=3.)
#    plot(29.98*bms.wns[bms.wsx],bms.ceffs,':',color='blue',linewidth=3.)
#    plot(29.98*bms.wnl[bms.wlx],bms.ceffl,':',color='blue',linewidth=3.)
#    plot(29.98*bms.wns[bms.wsx],bms.yfs*bms.ceffs,'--',color='black',linewidth=3.)
#    plot(29.98*bms.wnl[bms.wlx],bms.yfl*bms.ceffl,'--',color='black',linewidth=3.)
    
#    for i in range(4):
#        if i == 0:
#            pltfile='saturn_correct_v11_1.fits'
#            color='black'
#            label='Saturn'
#        if i == 1:
#            pltfile='mars_correct_v11_1.fits'
#            color='red'
#            label='Mars'
#        if i == 2:
#            pltfile='neptune_correct.fits'
#            color='blue'
#            label='Neptune'
#        if i == 3:
#            pltfile='uranus_correct.fits'
#            color='green'
#            label='Uranus'
#
#        pltfts=None
#        pltfts=pyfits.open(plotdir+pltfile)
#        temp=where(pltfts[2].data['flux'] > 1.2)
#        pltfts[2].data['flux'][temp]=1.
#        plot(pltfts[2].data['wave'],pltfts[2].data['flux'],color=color, \
#                 label=label)
#        temp=where(pltfts[4].data['flux'] > 1.2)
#        pltfts[4].data['flux'][temp]=1.
#        plot(pltfts[4].data['wave'],pltfts[4].data['flux'],color=color)
#        print label
#        print average(pltfts[2].data['flux']), \
#            average(pltfts[4].data['flux'])
#    errorbar(857.,0.69,xerr=0.,yerr=0.03,color='black')
#    errorbar(666.,0.59,xerr=0.,yerr=0.03,color='black')
#    plot(857.,0.69,marker='*',linewidth=0.,color='black')
#    plot(666.,0.59,marker='*',linewidth=0.,color='black')
    xlabel('Frequency (GHz)',fontsize='x-large')
    ylabel('$\\eta_{c}$',fontsize='x-large')
    xlim([400.,1600.])
    ylim([0.,0.8])
    legend(loc=1,fontsize='small',ncol=3)
    savefig(filename=plotdir+'gamma_all_v11.png',format='PNG')
    os.system('xview -shrink %sgamma_all_v11.png &' %plotdir)
    close(1)

def plot_planet(obj=None,caltype=None,nobs=None):

    sectdatapath='../data/FTS/sect_planets/'
    ftsdatapath='../data/FTS/'
    modelfile=None
    mc=False

    if obj == 'saturn':
        datapath=sectdatapath
        nobs=nobs
        d_plt_ini=16.
        obsid=['1342198279','1342224754','1342247750']
        dobs=[17.7,16.7,17.4]
        suffix=''
        modelfile='spire_JanskyGHz.dat'
        ftspsfile='./v11/Saturn/'+obsid[nobs]+ \
            '_spectrum_point_HR_unapod_18Apr2013' \
            +suffix+'.fits'
        dtrue=dobs[nobs]
        scale=(dtrue/dobs[0])**2.
        nskip=3
        munit=1.
        uplot=1e-3
        unit='kJy'
        yrange=[3.,40.]
        sim_func=sim_planet
        psuf='_v11_'+str(nobs)

    if obj == 'mars':
        datapath=sectdatapath
        d_plt_ini=6.
        nobs=nobs
        modelfile='MarsModel_OD176.txt'
        obsid=['1342193675', '1342197462', '1342198930', \
                   '1342231056', '1342231059', '1342231062', \
                   '1342231064', '1342231070', '1342231076', \
                   '1342231079', '1342231085', '1342245848', \
                   '1342247563', '1342247571']
        dobs=[8.8,6.0,5.4,5.5,5.5,5.5, \
                  5.5,5.5,5.5,5.5,5.5,8.8,6.6,6.6]
        ftspsfile='./v11/Mars/'+obsid[nobs]+'_spectrum_' \
            'point_HR_unapod_18Apr2013.fits'
        dtrue=dobs[nobs]
        scale=(dtrue/5.5)**2.
        nskip=2
        munit=1e3
        uplot=1e-3
        unit='kJy'
        yrange=[0.5,20.]
        sim_func=sim_planet
        psuf='_v11_'+str(nobs)

    if obj == 'neptune':
        datapath=sectdatapath
        d_plt_ini=2.
        dtrue=2.28
        modelfile='0x5000AD87_model_HR.fits'
        obsid='1342221703'
        ftspsfile='./v10/'+obsid+'_spectrum_' \
            'point_HR_unapod_v10.fits'
        scale=1.
        uplot=1.
        unit='Jy'
        yrange=[20.,300.]
        sim_func=sim_planet
        psuf=''

    if obj == 'uranus':
        datapath=sectdatapath
        d_plt_ini=3.
        dtrue=3.55
        modelfile='5001389B_model_HR.fits'
        obsid='1342257307'
        ftspsfile='./v10/'+obsid+'_spectrum_' \
            'point_HR_unapod_v10.fits'
        scale=1.
        uplot=1.
        unit='Jy'
        yrange=[90.,600.]
        sim_func=sim_planet
        psuf=''

    if obj == 'm83' or obj == 'M83':
        datapath=ftsdatapath
        modelfile=''
        d_plt_ini=18.
        ftspsfile='M83/HIPE.9.0.2742/unapodize/' \
            '1342212345_avgspectrum_HR_unapod_0_6.fits'
        ftsexfile='M83/HIPE.9.0.2742/unapodize/' \
            '1342212345_avgspectrum_HR_unapod_0_6.fits'
        yrange=[0.1,300.]
        uplot=1.
        unit='Jy'
        sim_func=sim_gauss
        psuf=''

    if obj == 'm82' or obj == 'M82':
        datapath=ftsdatapath
        modelfile=''
        d_plt_ini=22.
        ftspsfile='M82/HIPE.8.0.3384/unapodize/' \
            '1342208388_avgspectrum_HR_unapod_0_6.fits'
        ftsexfile='M82/HIPE.8.0.3384/unapodize/' \
            '1342208388_avgspectrum_HR_unapod_0_6.fits'
        yrange=[1.,2000.]
        uplot=1.
        unit='Jy'
        sim_func=sim_gauss
        psuf=''

    if obj == 'ngc4214' or obj == 'NGC4214':
        datapath=ftsdatapath
        modelfile=''
        d_plt_ini=18.
        ftspsfile=obj.upper()+'/HIPE.11.0.2785/unapodize/' \
            '1342256082_spectrum_point_HR_unapod.fits'
        ftsexfile='M83/HIPE.11.0.2785/unapodize/' \
            '1342256082_spectrum_extended_HR_unapod.fits'
        yrange=[0.01,20.]
        uplot=1.
        unit='Jy'
        sim_func=sim_gauss
        psuf=''

    if caltype == 'point':
        ftsfile=pyfits.open(datapath+ftspsfile)
    if caltype == 'extended':
        ftsfile=pyfits.open(datapath+ftsexfile)

    ftsfile.verify('silentfix')
    ctrxy=array([128,128])

    # ======= Read model spectrum of the planet =========

    if bool(modelfile):

        if modelfile[-4:] != 'fits':

            spec_mod=open('../data/FTS/sect_planets/'+modelfile)
            modline=spec_mod.readlines()
            smfreq=[]
            smkjy=[]

            for i in range(len(modline)/10-nskip):
                ii=i*10+nskip
                temp=None
                temp=modline[ii].split()
                smfreq.append(double(temp[0]))
                smkjy.append(double(temp[len(temp)-1]))
                mfreq=array(smfreq)
                mkjy=array(smkjy).copy()/munit

            if obj == 'mars':
                    
                mfreq=3e5/array(smfreq)
                scalem=0.
                
                for np in range(len(modelwn)):
                    
                    anp=abs(mfreq-modelwn[np]*30.)
                    idxnp=where(anp == min(anp))
                    scalem=scalem+(1e-3*(fluxLellouch[np]+fluxRudy[np])/2.)/ \
                        mkjy[idxnp[0][0]]
                    
                scalem=scalem/float(len(modelwn))
                
                mkjy=scalem*mkjy.copy()

        else:

            model=pyfits.open('../data/FTS/sect_planets/'+modelfile)
            model.verify('silentfix')
            nover=len((where(model['SSWD4'].data['wave'] <= \
                            model['SLWC3'].data['wave'].max()))[0])
            mslw=model['SLWC3'].data
            mssw=model['SSWD4'].data
            nmod=len(mslw['flux'])+len(mssw['flux'])-nover

            mfreq=zeros(nmod)
            mkjy=zeros(nmod)

            mfreq[0:len(mslw['flux'])]= \
                mslw['wave'][0:len(mslw['flux'])]
            mfreq[len(mslw['flux']):]=mssw['wave'][nover:]

            mkjy[0:(len(mslw['flux'])-nover)]= \
                mslw['flux'][0:(len(mslw['flux'])-nover)]
            mkjy[(len(mslw['flux'])-nover):len(mslw['flux'])]= \
                0.5*(mslw['flux'][(len(mslw['flux'])-nover):]+ \
                         mssw['flux'][0:nover])
            mkjy[len(mslw['flux']):]=mssw['flux'][nover:]
            

    # ======== Construct planet light distribution image =========

    [chiparam,chierr,dplt]=find_best_d(ftsfile,d_plt_ini, \
                                    caltype=caltype,obj=obj)
    prob=exp(-(chiparam-chiparam.min())/chiparam.min())
    mindx=where(chiparam == chiparam.min())
    mindx=mindx[0][0]
    proberr=zeros(len(prob))+1.
#    proberr=sqrt((chierr**2.+chiparam[mindx]**2.)) #\
#                     +(chiparam*chierr[mindx]/chiparam[mindx]**2.)**2.)

    figure(10)
    plot(dplt,prob,marker='*')
    fa={'x':dplt,'y':prob,'err':proberr}
    p_info=[{'value':d_plt_ini,'fixed':0, \
                 'limited':[1,1],'limits':[0.,100.]}, \
                {'value':2.,'fixed':0, \
                     'limited':[1,1],'limits':[0.,100.]}]
    params=mpfit.mpfit(fit_likelihood,functkw=fa,maxiter=1000,parinfo=p_info)
    dplt_g=dplt[0]+arange(101)*(dplt[-1]-dplt[0])/100.
    yy=fit_gaussian(dplt_g,params.params[0],params.params[1])
    plot(dplt_g,yy)
    xlim([dplt[0],dplt[-1]])
    ylim([0.,1.2])
    xlabel('Diameter (arcsec)')
    ylabel('exp(chi-param)')
    annotate('$\mu$ = {:.2f} (arcsec)'.format(params.params[0]), \
                 [dplt_g[0],1.1])
    annotate('$\sigma$ = {:.2f} (arcsec)'.format(params.params[1]), \
                 [dplt_g[0],1.0])
    savefig(plotdir+obj+'_d_fit'+psuf+'.png',format='PNG')
    close(10)

#    dplanet=params.params[0]
    dplanet=params.params[0]
    planet_mod=sim_func(dplanet,ctrxy)
    planet_area=sum(planet_mod)
    if obj == 'saturn':
        mod1=sim_gauss(12.,ctrxy)
        mod1_area=sum(mod1)
        mod2=sim_exp(4.,ctrxy)
        mod2_area=sum(mod2)

    # ======== Correcting the FTS Spectra ==========

    figure(0,figsize=[8.,8.])
    subplot(111)
    title('Correction on the '+obj.capitalize()+' FTS Spectrum',fontsize='large')
#    suptitle('$\\theta_{D}=$'+'{:.2f}'.format(params.params[0]) \
#                 +'   $\\theta_{obs}=$'+'{:.2f}'.format(dtrue))

    out=pyfits.HDUList()
    for i in range(2):
        if i == 0:
            wn=wn_ssw
            beam=beam_ssw
            fwhm=fwhm_ssw
            rtel=telrsrf[1]
            iro='blue'
            ni='SSW'
            deti='SSWD4'
        if i == 1:
            wn=wn_slw
            beam=beam_slw
            fwhm=fwhm_slw
            rtel=telrsrf[2]
            iro='red'
            ni='SLW'
            deti='SLWC3'

        fts=ftsfile[deti]

        if (obj == 'm83' or obj == 'M83') or \
                (obj == 'm82' or obj == 'M82'):
            fts.data['flux'][:]= \
                fts.data['flux'].copy()*cps[deti].data['pointConv']

        beam_area=[]
        beam_planet=[]
        if obj == 'saturn':
            beam_mod1=[]
            beam_mod2=[]
        for j in range(len(wn[0].data)):

            beam_area.append(sum(beam[0].data[j,:,:]))
            beam_planet.append(sum(beam[0].data[j,:,:]*planet_mod))
            if obj == 'saturn':
                beam_mod1.append(sum(beam[0].data[j,:,:]*mod1))
                beam_mod2.append(sum(beam[0].data[j,:,:]*mod2))

        a_beam=None
        ps_coup=None
        eta_mb=None
        a_beam=array(beam_area)

        if i == 0:
            assw=a_beam.copy()
        if i == 1:
            aslw=a_beam.copy()

        ps_coup=array(beam_planet)

        if obj == 'saturn':
            pmod1_coup=array(beam_mod1)
            pmod2_coup=array(beam_mod2)

        x=wn[0].data*30.   # change the wave number to GHz
        fab=interpolate.interp1d(x,a_beam, \
                                     bounds_error=False,fill_value=a_beam[0])
        fps=interpolate.interp1d(x,ps_coup, \
                                     bounds_error=False,fill_value=ps_coup[0])
        if obj == 'saturn':
            fmod1=interpolate.interp1d(x,pmod1_coup, \
                                         bounds_error=False,fill_value=pmod1_coup[0])
            fmod2=interpolate.interp1d(x,pmod2_coup, \
                                         bounds_error=False,fill_value=pmod2_coup[0])
            
        xp=ftsfile[deti].data['wave']
        ypa=ftsfile[deti].data['flux']*uplot
        rs=rtel.data['rsrf']/cps[deti].data['pointConv']* \
            fps(xp)/planet_area

        if caltype == 'point':
            correction=rtel.data['rsrf']/(rs*cps[deti].data['pointConv'])

        if caltype == 'extended':
            correction=rtel.data['rsrf']/rs

        ypb=ypa.copy()*correction
        if obj == 'saturn':
            ypb1=ypb/(planet_area/fps(xp)*fmod1(xp)/mod1_area)
            ypb2=ypb/(planet_area/fps(xp)*fmod2(xp)/mod2_area)
            fill_between(xp,ypb1,ypb2,color='gray',alpha=0.5)
            exec('foo'+str(i)+'=ypb1')
            exec('bar'+str(i)+'=ypb2')
        
        plot(xp,ypa,'--',color=iro,label=ni+' '+caltype,linewidth=1.)
        plot(xp,ypb,'-',color=iro,linewidth=1.,label=ni+' corrected')
        fts.data['flux'][:]=ypb[:]
        fts.data['error'][:]=fts.data['error']*correction
        out.append(fts)

        tempfts=fts.copy()

        if bool(modelfile):
            if modelfile[-4:] != 'fits':
                idx=(mfreq >=xp[0]-5.) & (mfreq <= xp[-1]+5.)
                fmod=interpolate.interp1d(mfreq[idx],mkjy[idx], \
                                              bounds_error=False, \
                                              fill_value=mkjy[-1])
                tempfts.data['flux'][:]=fmod(xp)*scale/ypb
            else:
                idx=(mfreq >= xp[0]) & (mfreq <= xp[-1])
                tempfts.data['flux'][:]=scale*mkjy[idx]/ypb

            tempfts.data['error'][:]=0.
            tempfts.data['mask'][:]=0.
            tempfts.name='GAMMA'+deti

            out.append(tempfts)


    if bool(modelfile):
        plot(mfreq,scale*mkjy,color='black',label='Model')
    if obj == 'mars':
        plot(modelwn*30.,scale*fluxLellouch*1e-3,linewidth=0.,marker='+', \
                 markersize=10.,color='black',markeredgewidth=2., \
                 label='Lellouch')
        plot(modelwn*30.,scale*fluxRudy*1e-3,linewidth=0.,marker='+', \
                 markersize=10.,color='red',markeredgewidth=2., \
                 label='Rudy')
    
    xlabel('Frequency (GHz)',fontsize='x-large')
    ylabel('Flux Density('+unit+')',fontsize='x-large')
    xlim([400.,1600.])
    ylim(yrange)
#    annotate('{:.2f}'.format(params.params[0]),[450.,0.95*yrange[0]])
    axvline(linewidth=3.,color='black')
    axhline(linewidth=3.,color='black')
    yscale('log')
    legend(loc=2)
    savefig(filename=plotdir+obj+psuf+'.png',format='PNG')
    close(0)    

    os.system('open %s%s%s.png &'%(plotdir,obj,psuf))
    out.writeto(plotdir+obj+'_correct'+psuf+'.fits',clobber=True, \
                    output_verify='ignore')

    return

def cal_fbeamarea(i):

    asc2sr=(pi/180./3600.)**2.

    if i == 'ssw':
        beam=beam_ssw
        wn=wn_ssw[0].data[:]*30.

    if i == 'slw':
        beam=beam_slw
        wn=wn_slw[0].data[:]*30.

    beam_area=[]

    for j in range(len(wn)):
        beam_area.append(sum(beam[0].data[j,:,:]))

    a_beam=None
    a_beam=array(beam_area)

    if i == 'ssw':
        area=a_beam.copy()*asc2sr
        f_a=interpolate.interp1d(wn,area,kind='cubic',fill_value=area[0], \
                                           bounds_error=False)
    if i == 'slw':
        area=a_beam.copy()*asc2sr
        f_a=interpolate.interp1d(wn,area,kind='cubic',fill_value=area[0], \
                                           bounds_error=False)

    return f_a
 
def cs_vs_size():

    ssw_ov=where(wn_ssw[0].data*30. < fmax_slw)
    slw_ov=where(wn_slw[0].data*30. > fmin_ssw)
    ssw_ov=ssw_ov[0][:]
    slw_ov=slw_ov[0][0:-1]

    ctrxy=[128.,128.]
    dsim=arange(61)*2.

    for pp in range(3):

        if pp == 0:
            model=sim_gauss
            rin=dsim.copy()
        if pp == 1:
            model=sim_exp
            rin=dsim.copy()
        if pp == 2:
            model=sim_planet
            rin=dsim.copy()

        output=[]
        output.append(0.)
        for i in range(len(rin)-1):
            profile=model(rin[i+1],ctrxy)

            fac_ssw=[]
            fac_slw=[]

            for j in range(len(ssw_ov)):
                fac_ssw.append(sum(profile)/ \
                                   sum(beam_ssw[0].data[ssw_ov[j],:,:]*profile))
                fac_slw.append(sum(profile)/ \
                                   sum(beam_slw[0].data[slw_ov[j],:,:]*profile))

            fac_ssw=array(fac_ssw)
            fac_slw=array(fac_slw)

            output.append(average((fac_ssw-fac_slw)/fac_ssw))

        if pp == 0:
            output_gp=array(output)
        if pp == 1:
            output_ep=array(output)
        if pp == 2:
            output_pp=array(output)

    maxvalue=average(1- \
                         sum(sum(beam_ssw[0].data[ssw_ov,:,:],axis=1),axis=1)/ \
                         sum(sum(beam_slw[0].data[slw_ov,:,:],axis=1),axis=1))
    y1=1.-(1.-maxvalue)*1.197718322035786
    y2=1.-(1.-maxvalue)*0.6822450815009493
    sw1d=pyfits.open('../data/FTS/FTS_Beam/SSW_beam_1d.2.fits')
    lw1d=pyfits.open('../data/FTS/FTS_Beam/SLW_beam_1d.2.fits')
    bmax=lw1d[1].data['fwhm'].max()
    bmin=sw1d[1].data['fwhm'].min()

    figure(0,figsize=[8.,6.])
#    title('$\\Delta F/F_{SLW}$ at overlap')
    plot(dsim,(zeros(len(dsim))+1.)*maxvalue*100.,'--',color='black')
    fill_between(dsim,(zeros(len(dsim))+1.)*y1*100., \
                     (zeros(len(dsim))+1.)*y2*100., \
                     color='grey',alpha=0.2)
    plot(bmax*(zeros(2)+1.),[0.,100.],'--',color='black')
    plot(bmin*(zeros(2)+1.),[0.,100.],'--',color='black')
    print bmax, bmin
    plot(dsim,output_gp*100.,linewidth=3.,label='Gaussian',color='red')
    plot(dsim,output_ep*100.,linewidth=3.,label='Sersic, n=1',color='blue')
    plot(dsim,output_pp*100.,linewidth=3.,label='Top Hat',color='green')
    xlabel('$\\theta_{D}$ (")',fontsize='x-large')
    xlim([0.,120.])
    ylabel('$\\frac{(F_{SLW}-F_{SSW})}{F_{SLW}}$ (%)',fontsize='x-large')
    ylim([0.,100.])
    legend(loc=4)
    savefig(filename='../results/SECT/overlap_size.png',format='PNG')
    close(0)

def maxb_mask(dfwhm):

    rcut=3.*dfwhm/(2.*sqrt(2.*log(2.)))
    img_mask=zeros([257,257])+1.
    imgsize=img_mask.shape
    mask_ctr=int(round(imgsize[0]*0.5))
    x=arange(imgsize[0])
    y=arange(imgsize[1])
    for xi in range(mask_ctr,imgsize[0]):
        idx=where(((x[xi]-mask_ctr)**2.+(y-mask_ctr)**2.) > (rcut)**2.)
        img_mask[x[xi],y[idx]]=0.
        img_mask[2.*mask_ctr-x[xi],y[idx]]=0.

    return img_mask

# ======== Common block to load ========

side_det=[['SSWB2','SSWB3','SSWB4', \
               'SSWC2','SSWC3','SSWC4','SSWC5', \
               'SSWD2','SSWD3','SSWD4','SSWD6', \
               'SSWE2','SSWE3','SSWE4','SSWE5', \
               'SSWF2','SSWF3'], \
              ['SLWB2','SLWB3', \
                   'SLWC2','SLWC3','SLWC4', \
                   'SLWD2','SLWD3']]

fts_beamdir='../data/FTS/FTS_Beam/'
caldir='../data/FTS/spire_cal_12_3/'
plotdir='../results/SECT/'
wn_slw=pyfits.open(fts_beamdir+'slw_wn.2.fits')
wn_slw.verify('silentfix')
wn_ssw=pyfits.open(fts_beamdir+'ssw_wn.2.fits')
wn_ssw.verify('silentfix')
beam_slw=pyfits.open(fts_beamdir+'slw_beam.2.fits')
beam_slw.verify('silentfix')
beam_ssw=pyfits.open(fts_beamdir+'ssw_beam.2.fits')
beam_ssw.verify('silentfix')
fwhm_slw=pyfits.open(fts_beamdir+'SLW_beam_1d.2.fits')
fwhm_slw.verify('silentfix')
fwhm_ssw=pyfits.open(fts_beamdir+'SSW_beam_1d.2.fits')
fwhm_ssw.verify('silentfix')
cps=pyfits.open(caldir+'herschel.spire.ia.dataset.SpecBeamParam/' \
                    'SCalSpecBeamParam_HR_20120218_v3.fits')
cps.verify('silentfix')
rtel=pyfits.open(caldir+'herschel.spire.ia.dataset.SpecTeleRsrf/' \
                     'SCalSpecTeleRsrf_6_HR_20120218_v2.fits')
rtel.verify('silentfix')

img_mask=maxb_mask(42.)

if os.path.exists(plotdir+'telersrf_avg.fits'):
    telrsrf=pyfits.open(plotdir+'telersrf_avg.fits')
    telrsrf.verify('silentfix')
if os.path.exists(plotdir+'bol_vresp.fits'):
    bol_vres=pyfits.open(plotdir+'bol_vresp.fits')
    bol_vres.verify('silentfix')
else:
    save_rsrf_fits()

fmax_slw=cps['SLWC3'].data['wave'].max()
fmin_ssw=cps['SSWD4'].data['wave'].min()

# Mars continuum model from Sunil for 0x50005492
modelwn   = array([20, 30, 40, 50])
fluxLellouch = array([1145.8, 2520.4, 4361.1, 6617.4])
fluxRudy     = array([1151.54, 2510.79, 4322.45, 6538.06])
 
# Mars continuum model from Sunil for 0x50004ED6
#fluxLellouch = array([1395.50, 3067.58, 5305.91, 8048.52])                 
#fluxRudy     = array([1410.51, 3074.72, 5292.04, 8002.71])


# correct errors in fitting

# flux=
# errors=array([  88.72194807,  117.14844707,  103.30168961,  125.27772113,
         #94.39458379,  120.01647721,   73.50222375,   56.06274772,
         #40.25470352,  123.68647411])

# wrong errors in fitting

# flux=

# errors=
