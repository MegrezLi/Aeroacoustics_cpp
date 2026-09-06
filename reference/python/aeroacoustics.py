# Copyright (C) 2012-2016 National Renewable Energy Laboratory
# Derived from OpenFAST (Apache-2.0). Python modifications: 2026-09-05.
# Licensed under the Apache License, Version 2.0. See LICENSE-OpenFAST.txt.
# Provided AS IS, without warranties or conditions of any kind.
"""OpenFAST aeroacoustics in Python/NumPy; no executable, DLL or Fortran runtime.

Source: OpenFAST commit 2895884d2be01862173c88d70f86b358d2f1a50a.
Run: python aeroacoustics.py --output results
All kernel names correspond to their Fortran counterparts (lowercase).
Units: SI, SPL in dB re 20 microPa, all angles degrees except inflownoise alpha.
The default section example adapts the numeric case in AeroAcoustics.f90:Aero_Tests.
This is the acoustic module, not the structural/aerodynamic OpenFAST solver.
"""
from __future__ import annotations

import argparse
import csv
import json
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path
import numpy as np

SOURCE_COMMIT = '2895884d2be01862173c88d70f86b358d2f1a50a'
aa_epsilon = 1e-16
x_blmethod_tables = 2
itrip_none, itrip_heavy, itrip_light = 0, 1, 2
pi, twopi = np.pi, 2 * np.pi
FREQUENCIES = np.array([10,12.5,16,20,25,31.5,40,50,63,80,
    100,125,160,200,250,315,400,500,630,800,1000,1250,1600,2000,
    2500,3150,4000,5000,6300,8000,10000,12500,16000,20000], dtype=float)
MECHANISMS = ('laminar', 'pressure', 'suction', 'separation', 'blunt', 'tip', 'inflow')


@dataclass
class Parameters:
    """Acoustic parameters; switches follow the OpenFAST input file numbering."""
    freqlist: np.ndarray = field(default_factory=lambda: FREQUENCIES.copy())
    spdsound: float = 340.0
    kinvisc: float = 1.48e-5
    airdens: float = 1.225
    lturb: float = 40.0
    x_blmethod: int = 1
    itrip: int = 1
    round: bool = True
    alprat: float = 1.0
    timod: int = 1
    tbltemod: int = 1
    lammod: int = 0
    tipmod: int = 0
    bluntmod: int = 0
    aweighting: bool = False
    ti: float = 0.1
    avgv: float = 8.0

    def __post_init__(self):
        self.freqlist = np.asarray(self.freqlist, dtype=float)
        if self.freqlist.ndim != 1 or len(self.freqlist) == 0 or not np.all(np.isfinite(self.freqlist)) or np.any(self.freqlist <= 0) or np.any(np.diff(self.freqlist) <= 0):
            raise ValueError('freqlist must be finite, positive and strictly increasing')
        for name in ('spdsound', 'kinvisc', 'airdens', 'lturb'):
            if not np.isfinite(getattr(self, name)) or getattr(self, name) <= 0:
                raise ValueError(f'{name} must be finite and positive')
        for name, allowed in {'timod': (0,1,2), 'tbltemod': (0,1,2),
                              'x_blmethod': (1,2), 'itrip': (0,1,2),
                              'lammod': (0,1), 'tipmod': (0,1), 'bluntmod': (0,1)}.items():
            if getattr(self,name) not in allowed:
                raise ValueError(f'{name} must be one of {allowed}')
        if not np.all(np.isfinite([self.ti,self.avgv,self.alprat])) or self.ti < 0 or self.avgv < 0 or self.alprat <= 0:
            raise ValueError('ti/avgv must be nonnegative; alprat must be positive')


def a_weighting(frequency):
    """Exact SetParameters expression, including upstream rounded constants."""
    f2 = np.asarray(frequency, dtype=float)**2
    f4 = f2*f2
    return (10*np.log10(1.562339*f4/((f2+107.65265**2)*(f2+737.86223**2)))
            + 10*np.log10(2.242881e16*f4/((f2+20.598997**2)**2*(f2+12194.22**2)**2)))


def db_sum(levels, axis=None):
    """Sum incoherent sources in energy space. Silence is -inf, never 0 dB."""
    levels = np.asarray(levels, dtype=float)
    with np.errstate(divide='ignore'):
        return 10/np.log(10)*np.logaddexp.reduce(levels*np.log(10)/10, axis=axis)

# BEGIN PORTED FORTRAN KERNELS

def log10aa(x):
    """Port of LOG10AA, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns f."""
    f = np.log10( max(aa_epsilon, x) )
    return f

def lblvs(alpstar, c, u, theta, phi, l, r, p, d99var2, dstarvar1, dstarvar2, stallval):
    """Port of LBLVS, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns spllam."""
    spllam = np.zeros(len(p.freqlist))
    m = u  / p.spdsound
    rc = u  * c/p.kinvisc
    if p.x_blmethod  ==  x_blmethod_tables:
        deltap = d99var2
        dstrs = dstarvar1
        dstrp = dstarvar2
    else:
        deltap, dstrs, dstrp = thick(c, rc, alpstar, p, stallval)
    dbarh = directh_te(m,theta,phi)
    if dbarh <= 0:
        spllam[:] = 0.
        return spllam
    if rc  <=  1.3e+05:
        st1prim = .18
    elif rc <= 4.0e+05:
        st1prim = .001756*rc**.3931
    else:
        st1prim = .28
    stpkprm = 10.**(-.04*alpstar) * st1prim
    if alpstar  <=  3.0:
        rc0 = 10.**(.215*alpstar+4.978)
    else:
        rc0 = 10.**(.120*alpstar+5.263)
    d = rc / rc0
    if d  <=  .3237:
        g2 = 77.852*log10aa(d)+15.328
    elif d  <=  .5689:
        g2 = 65.188*np.log10(d) + 9.125
    elif d  <=  1.7579:
        g2 = -114.052 * np.log10(d)**2
    elif d  <=  3.0889:
        g2 = -65.188*np.log10(d)+9.125
    else:
        g2 = -77.852*np.log10(d)+15.328
    g3 = 171.04 - 3.03 * alpstar
    scale = 10. * log10aa(deltap*m**5*dbarh*l/r**2)
    for i in range(1, len(p.freqlist) + 1):
        stprim = p.freqlist[i - 1] * deltap / u
        e = stprim / stpkprm
        if e  <=  .5974:
            g1 = 39.8*log10aa(e)-11.12
        elif e  <=  .8545:
            g1 = 98.409 * np.log10(e) + 2.0
        elif e  <=  1.17:
            g1 = -5.076+np.sqrt(2.484-506.25*(np.log10(e))**2)
        elif e  <=  1.674:
            g1 = -98.409 * np.log10(e) + 2.0
        else:
            g1 = -39.80*np.log10(e)-11.12
        spllam[i - 1] = g1 + g2 + g3 + scale
    return spllam

def tblte(alpstar, c, u, theta, phi, l, r, p, d99var2, dstarvar1, dstarvar2, stallval):
    """Port of TBLTE, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns splp, spls, splalph."""
    splp = np.zeros(len(p.freqlist))
    spls = np.zeros(len(p.freqlist))
    splalph = np.zeros(len(p.freqlist))
    m = u  / p.spdsound
    rc = u  * c/p.kinvisc
    if p.x_blmethod  ==  x_blmethod_tables:
        deltap = d99var2
        dstrs = dstarvar1
        dstrp = dstarvar2
    else:
        deltap, dstrs, dstrp = thick(c, rc, alpstar, p, stallval)
    dbarl = directl(m,theta,phi)
    dbarh = directh_te(m,theta,phi)
    rdstrs = dstrs * u  / p.kinvisc
    rdstrp = dstrp * u  / p.kinvisc
    st1 = .02 * m ** (-.6)
    if alpstar  <=  1.333:
        st2 = st1
    elif alpstar  <=  stallval:
        st2 = st1*10.**(.0054*(alpstar-1.333)**2)
    else:
        st2 = 4.72 * st1
    st1prim = (st1+st2)/2.
    a0 = a0comp(rc)
    a02 = a0comp(3.0*rc)
    amina0 = amin(a0)
    amaxa0 = amax(a0)
    amina02 = amin(a02)
    amaxa02 = amax(a02)
    ara0 = (20. + amina0) / (amina0 - amaxa0)
    ara02 = (20. + amina02)/ (amina02- amaxa02)
    if rc  <  9.52e+04:
        b0 = .30
    elif rc  <  8.57e+05:
        b0 = (-4.48e-13)*(rc-8.57e+05)**2 + .56
    else:
        b0 = .56
    bminb0 = bmin(b0)
    bmaxb0 = bmax(b0)
    brb0 = (20. + bminb0) / (bminb0 - bmaxb0)
    stpeak = st1
    if rc  <  2.47e+05:
        k1 = -4.31 * log10aa(rc) + 156.3
    elif rc  <=  8.0e+05:
        k1 = -9.0 * np.log10(rc) + 181.6
    else:
        k1 = 128.5
    if rdstrp  <=  5000.:
        delk1 = -alpstar*(5.29-1.43*log10aa(rdstrp))
    else:
        delk1 = 0.0
    gamma = 27.094 * m +  3.31
    beta = 72.650 * m + 10.74
    gamma0 = 23.430 * m +  4.651
    beta0 = -34.190 * m - 13.820
    if alpstar  <=  (gamma0-gamma):
        k2 = -1000.0
    elif alpstar <= (gamma0+gamma):
        k2 = np.sqrt(beta**2-(beta/gamma)**2*(alpstar-gamma0)**2)+beta0
    else:
        k2 = -12.0
    k2 = k2 + k1
    xcheck = gamma0
    switch = False
    if (alpstar  >=  xcheck) or (alpstar  >  stallval):
        switch = True
    for i in range(1, len(p.freqlist) + 1):
        stp = p.freqlist[i - 1] * dstrp / u
        a = log10aa( stp / stpeak )
        amina = amin(a)
        amaxa = amax(a)
        aa = amina + ara0 * (amaxa - amina)
        splp[i - 1] = aa+k1-3.+10.*log10aa(dstrp*m**5*dbarh*l/r**2)+delk1
        sts = p.freqlist[i - 1] * dstrs / u
        if not  switch:
            a = log10aa( sts / st1prim )
            amina = amin(a)
            amaxa = amax(a)
            aa = amina + ara0 * (amaxa - amina)
            spls[i - 1] = aa+k1-3.+10.*log10aa(dstrs*m**5*dbarh* l/r**2)
            b = log10aa(sts / st2)
            bminb = bmin(b)
            bmaxb = bmax(b)
            bb = bminb + brb0 * (bmaxb-bminb)
            splalph[i - 1] = bb+k2+10.*log10aa(dstrs*m**5*dbarh*l/r**2)
        else:
            spls[i - 1] = 10.*log10aa(dstrs*m**5*dbarl*l/r**2)
            splp[i - 1] = 10.*log10aa(dstrp*m**5*dbarl*l/r**2)
            b = log10aa(sts / st2)
            aminb = amin(b)
            amaxb = amax(b)
            bb = aminb + ara02 * (amaxb-aminb)
            splalph[i - 1] = bb+k2+10.*log10aa(dstrs*m**5*dbarl*l/r**2)
        if splp[i - 1]     <  -100.:
            splp[i - 1] = -100.
        if spls[i - 1]     <  -100.:
            spls[i - 1] = -100.
        if splalph[i - 1]  <  -100.:
            splalph[i - 1] = -100.
    return splp, spls, splalph

def tipnois(alphtip, alprat2, c, u, theta, phi, r, p):
    """Port of TIPNOIS, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns spltip."""
    spltip = np.zeros(len(p.freqlist))
    if alphtip == 0.:
        spltip[:] = 0
        return spltip
    elif alphtip < 0.:
        pass
    alptipp = abs(alphtip) * alprat2
    m = u  / p.spdsound
    dbarh = directh_te(m,theta,phi)
    if p.round:
        l = .008 * alptipp * c
    else:
        if abs(alptipp)  <=  2.:
            l = (.023 + .0169*alptipp) * c
        else:
            l = (.0378 + .0095*alptipp) * c
    mm = (1. + .036*alptipp) * m
    um = mm * p.spdsound
    term = m*m*mm**3*l**2*dbarh/r**2
    if term  !=  0.0:
        scale = 10.*np.log10(term)
    else:
        scale = 0.0
    for i in range(1, len(p.freqlist) + 1):
        stpp = p.freqlist[i - 1] * l / um
        spltip[i - 1] = 126.-30.5*(log10aa(stpp)+.3)**2 + scale
    return spltip

def inflownoise(alphanoise, chord, u, theta, phi, d, robs, tinoise, p):
    """Port of INFLOWNOISE, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns splti."""
    splti = np.zeros(len(p.freqlist))
    mach = u/p.spdsound
    tinooisess = tinoise
    dbarl = directl(mach,theta,phi)
    dbarh = directh_le(mach,theta,phi)
    if dbarh <= 0:
        splti[:] = 0.
        return splti
    frequency_cutoff = 10*u/pi/chord
    ke = 3.0/(4.0*p.lturb)
    beta2 = 1-mach*mach
    for i in range(1, len(p.freqlist) + 1):
        if p.freqlist[i - 1] <= frequency_cutoff:
            directivity = dbarl
        else:
            directivity = dbarh
        wavenumber = twopi*p.freqlist[i - 1]/u
        kbar = wavenumber*chord/2.0
        khat = wavenumber/ke
        splhigh = 10.*log10aa(p.airdens**2 * p.spdsound**4 * p.lturb * (d/2.) / (robs**2) *(mach**5) * tinooisess**2 *(khat**3)* (1+khat**2)**(-7./3.) * directivity) + 78.4
        splhigh = splhigh + 10.*np.log10(1+ 9.0*alphanoise**2)
        sears = 1./(twopi*kbar/beta2 + 1./(1.+2.4*kbar/beta2))
        lfc = max(aa_epsilon, 10*sears*mach*kbar**2/beta2)
        splti[i - 1] = splhigh + 10.*log10aa(lfc/(1+lfc))
    return splti

def blunt(alpstar, c, u, theta, phi, l, r, h, psi, p, d99var2, dstarvar1, dstarvar2, stallval):
    """Port of BLUNT, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns splblunt."""
    splblunt = np.zeros(len(p.freqlist))
    g5 = np.zeros(len(p.freqlist))
    m = u  / p.spdsound
    rc = u  * c/p.kinvisc
    if p.x_blmethod  ==  x_blmethod_tables:
        deltap = d99var2
        dstrs = dstarvar1
        dstrp = dstarvar2
    else:
        deltap, dstrs, dstrp = thick(c, rc, alpstar, p, stallval)
    dstravg = (dstrs + dstrp) / 2.
    hdstar = h / dstravg
    dstarh = 1. /hdstar
    dbarh = directh_te(m,theta,phi)
    if dbarh <= 0:
        splblunt[:] = 0.
        return splblunt
    aterm = .212 - .0045 * psi
    if hdstar  >=  .2:
        stpeak = aterm / (1.+.235*dstarh-.0132*dstarh**2)
    else:
        stpeak = .1 * hdstar + .095 - .00243 * psi
    if hdstar  <=  5.:
        g4 = 17.5*log10aa(hdstar)+157.5-1.114*psi
    else:
        g4 = 169.7 - 1.114 * psi
    scale = 10. * log10aa(m**5.5 * h * dbarh * l / r**2)
    g5sum = 0.0
    for i in range(1, len(p.freqlist) + 1):
        stppp = p.freqlist[i - 1] * h / u
        eta = log10aa(stppp/stpeak)
        g514 = g5comp(hdstar,eta)
        hdstarp = 6.724 * hdstar **2-4.019*hdstar+1.107
        g50 = g5comp(hdstarp,eta)
        g5[i - 1] = g50 + .0714 * psi * (g514-g50)
        if g5[i - 1]  >  0.:
            g5[i - 1] = 0.
        g5sum = 10**(g5[i - 1]/10)+g5sum
        if g5sum  !=  0:
            logval = max(aa_epsilon,1/g5sum)
        else:
            logval = 1
        splblunt[i - 1] = g4 + g5[i - 1] + scale - 10*np.log10(logval)
    return splblunt

def g5comp(hdstar, eta):
    """Port of G5COMP, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns g5."""
    if hdstar  <  .25:
        mu = .1211
    elif hdstar  <=  .62:
        mu = -.2175*hdstar + .1755
    elif hdstar  <  1.15:
        mu = -.0308*hdstar + .0596
    else:
        mu = .0242
    if hdstar  <=  .02:
        m = 0.0
    elif hdstar  <  0.5:
        m = 68.724*hdstar - 1.35
    elif hdstar  <=  .62:
        m = 308.475*hdstar - 121.23
    elif hdstar  <=  1.15:
        m = 224.811*hdstar - 69.354
    elif hdstar  <  1.2:
        m = 1583.28*hdstar - 1631.592
    else:
        m = 268.344
    m = max(m, 0.0)
    eta0 = -np.sqrt((m*m*mu**4)/(6.25+m*m*mu*mu))
    if eta  <=  eta0:
        k = 2.5*np.sqrt(1.-(eta0/mu)**2)-2.5-m*eta0
        g5 = m * eta + k
    elif eta  <=  0.:
        g5 = 2.5*np.sqrt(1.-(eta/mu)**2)-2.5
    elif eta  <=  0.03615995:
        g5 = np.sqrt(1.5625-1194.99*eta**2)-1.25
    else:
        g5 = -155.543 * eta + 4.375
    return g5

def amin(a):
    """Port of AMIN, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns amina."""
    x1 = abs(a)
    if x1  <=  .204:
        amina = np.sqrt(67.552-886.788*x1**2)-8.219
    elif x1  <=  .244:
        amina = -32.665*x1+3.981
    else:
        amina = -142.795*x1**3+103.656*x1**2-57.757*x1+6.006
    return amina

def amax(a):
    """Port of AMAX, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns amaxa."""
    x1 = abs(a)
    if x1  <=  .13:
        amaxa = np.sqrt(67.552-886.788*x1**2)-8.219
    elif x1  <=  .321:
        amaxa = -15.901*x1+1.098
    else:
        amaxa = -4.669*x1**3+3.491*x1**2-16.699*x1+1.149
    return amaxa

def bmin(b):
    """Port of BMIN, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns bminb."""
    x1 = abs(b)
    if x1  <=  .13:
        bminb = np.sqrt(16.888-886.788*x1**2)-4.109
    elif x1  <=  .145:
        bminb = -83.607*x1+8.138
    else:
        bminb = -817.81*x1**3+355.21*x1**2-135.024*x1+10.619
    return bminb

def bmax(b):
    """Port of BMAX, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns bmaxb."""
    x1 = abs(b)
    if x1  <=  .1:
        bmaxb = np.sqrt(16.888-886.788*x1**2)-4.109
    elif x1  <=  .187:
        bmaxb = -31.313*x1+1.854
    else:
        bmaxb = -80.541*x1**3+44.174*x1**2-39.381*x1+2.344
    return bmaxb

def a0comp(rc):
    """Port of A0COMP, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns a0."""
    if rc  <  9.52e+04:
        a0 = .57
    elif rc  <  8.57e+05:
        a0 = (-9.57e-13)*(rc-8.57e+05)**2 + 1.13
    else:
        a0 = 1.13
    return a0

def thick(c, rc, alpstar, p, stallval):
    """Port of THICK, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns deltap, dstrs, dstrp."""
    logrc = log10aa( rc )
    delta0 = 10.**(1.6569-0.9045*logrc+0.0596*logrc**2)*c
    if p.itrip != itrip_none:
        delta0 = 10.**(1.892 -0.9045*logrc+0.0596*logrc**2)*c
    if p.itrip  ==  itrip_light:
        delta0 = .6*delta0
    deltap = 10.**(-.04175*alpstar+.00106*alpstar**2)*delta0
    if p.itrip != itrip_none:
        if rc  <=  .3e+06:
            dstr0 = .0601 * rc **(-.114)*c
        else:
            dstr0 = 10.**(3.411-1.5397*logrc+.1059*logrc**2)*c
        if p.itrip  ==  itrip_light:
            dstr0 = dstr0 * .6
    else:
        dstr0 = 10.**(3.0187-1.5397*logrc+.1059*logrc**2)*c
    dstrp = 10.**(-.0432*alpstar+.00113*alpstar**2)*dstr0
    if p.itrip  ==  itrip_heavy:
        if alpstar  <=  5.:
            dstrs = 10.**(.0679*alpstar)*dstr0
        elif alpstar  <=  stallval:
            dstrs = 0.381 * 10.**(.1516*alpstar)*dstr0
        else:
            dstrs = 14.296 * 10.**(.0258*alpstar)*dstr0
    else:
        if alpstar  <=  7.5:
            dstrs = 10.**(.0679*alpstar)*dstr0
        elif alpstar  <=  stallval:
            dstrs = .0162*10.**(.3066*alpstar)*dstr0
        else:
            dstrs = 52.42*10.**(.0258*alpstar)*dstr0
    return deltap, dstrs, dstrp

def directh_te(m, theta, phi):
    """Port of DIRECTH_TE, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns dbar."""
    degrad = .017453
    mc = .8 * m
    thetar = theta * degrad
    phir = phi * degrad
    dbar = 2.*np.sin(thetar/2.)**2 * np.sin(phir)**2 / ((1.+m*np.cos(thetar))* (1.+(m-mc)*np.cos(thetar))**2)
    return dbar

def directh_le(m, theta, phi):
    """Port of DIRECTH_LE, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns dbar."""
    degrad = .017453
    thetar = theta * degrad
    phir = phi * degrad
    dbar = 2.*np.cos(thetar/2.)**2*np.sin(phir)**2/(1.+m*np.cos(thetar))**3
    return dbar

def directl(m, theta, phi):
    """Port of DIRECTL, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns dbar."""
    degrad  = .017453
    mc = .8 * m
    thetar = theta * degrad
    phir = phi * degrad
    dbar = (np.sin(thetar)*np.sin(phir))**2/(1.+m*np.cos(thetar))**4
    return dbar

def simple_guidati(u, chord, thick_10p, thick_1p, p):
    """Port of SIMPLE_GUIDATI, AeroAcoustics.f90; SI units, angles in degrees
    (inflow_noise/Inflownoise alpha is radians). Returns splti."""
    splti = np.zeros(len(p.freqlist))
    ti_param = thick_1p + thick_10p
    slope = 1.123*ti_param + 5.317*ti_param*ti_param
    const1 = -slope*twopi*chord/u
    const2 = -slope*5.0e0
    for loop1 in range(1, len(p.freqlist) + 1):
        splti[loop1 - 1] = const1 * p.freqlist[loop1 - 1] + const2
    return splti


# Kronrod nodes/weights from SLATEC DQK61 (Piessens & de Doncker).
_XGK = np.array([0.9994844100504906, 0.9968934840746495, 0.9916309968704046, 0.9836681232797472, 0.9731163225011262, 0.9600218649683075, 0.94437444474856, 0.9262000474292743, 0.9055733076999078, 0.8825605357920527, 0.8572052335460612, 0.8295657623827684, 0.799727835821839, 0.7677774321048262, 0.7337900624532268, 0.6978504947933158, 0.6600610641266269, 0.6205261829892429, 0.5793452358263617, 0.5366241481420199, 0.49248046786177857, 0.44703376953808915, 0.4004012548303944, 0.3527047255308781, 0.30407320227362505, 0.25463692616788985, 0.20452511668230988, 0.15386991360858354, 0.10280693796673702, 0.0514718425553177, 0.0])
_WGK = np.array([0.0013890136986770077, 0.003890461127099884, 0.0066307039159312926, 0.009273279659517764, 0.011823015253496341, 0.014369729507045804, 0.01692088918905327, 0.019414141193942382, 0.021828035821609193, 0.0241911620780806, 0.0265099548823331, 0.02875404876504129, 0.030907257562387762, 0.03298144705748372, 0.034979338028060025, 0.03688236465182123, 0.038678945624727595, 0.040374538951535956, 0.041969810215164244, 0.04345253970135607, 0.04481480013316266, 0.04605923827100699, 0.04718554656929915, 0.04818586175708713, 0.04905543455502978, 0.04979568342707421, 0.05040592140278235, 0.05088179589874961, 0.051221547849258774, 0.05142612853745902, 0.05149472942945157])
_WG = np.array([0.007968192496166605, 0.01846646831109096, 0.02878470788332337, 0.03879919256962705, 0.04840267283059405, 0.057493156217619065, 0.06597422988218049, 0.0737559747377052, 0.08075589522942021, 0.08689978720108298, 0.09212252223778612, 0.09636873717464425, 0.09959342058679527, 0.1017623897484055, 0.10285265289355884])


def qk61(function, lower, upper, diagnostics=False):
    """SLATEC fixed 61-point Gauss-Kronrod rule in Python.

    A single panel, as used by OpenFAST; this is not adaptive quadrature.
    function accepts a NumPy vector. diagnostics returns value/error/abs/asc.
    """
    centre, half = (lower+upper)/2, (upper-lower)/2
    fc = function(np.asarray(centre))
    fm, fp = function(centre-half*_XGK[:-1]), function(centre+half*_XGK[:-1])
    resk = _WGK[-1]*fc + np.sum(_WGK[:-1]*(fm+fp))
    result = resk*half
    if not diagnostics:
        return result
    resg = np.sum(_WG*(fm[1::2]+fp[1::2]))
    resabs = (abs(fc)*_WGK[-1]+np.sum(_WGK[:-1]*(abs(fm)+abs(fp))))*abs(half)
    mean = resk/2
    resasc = (_WGK[-1]*abs(fc-mean)+np.sum(_WGK[:-1]*(abs(fm-mean)+abs(fp-mean))))*abs(half)
    error = abs((resk-resg)*half)
    if resasc != 0 and error != 0:
        error = resasc*min(1., (200*error/resasc)**1.5)
    eps, tiny = np.finfo(float).eps, np.finfo(float).tiny
    if resabs > tiny/(50*eps):
        error = max(50*eps*resabs, error)
    return result, error, resabs, resasc


def spl_integrate(omega, limits, issuctionside, mach, spdsound, airdens,
                  kinvisc, cfall, d99all, edgevelall):
    """TNO:SPL_integrate, Pressure, f_int1/f_int2; no shared global state.

    Two nested QK61 panels. Both sides are ordered [suction, pressure].
    The outer and inner evaluations are broadcast to a 61x61 NumPy array.
    """
    side = 0 if issuctionside else 1
    cf, delta = cfall[side], d99all[side]
    if cf <= 0 or delta <= 0:
        raise ValueError('TNO integral requires positive Cf and boundary-layer thickness')
    uo = mach*spdsound*abs(edgevelall[side])
    if uo <= 0 or kinvisc <= 0:
        raise ValueError('TNO integral requires positive edge speed and viscosity')
    # Explicit node set includes the midpoint exactly once.
    nodes = np.concatenate((-_XGK[:-1], [0.], _XGK[:-1]))
    weights = np.concatenate((_WGK[:-1], _WGK[-1:], _WGK[:-1]))
    k1 = ((limits[0]+limits[1])/2+(limits[1]-limits[0])/2*nodes)[:, None]
    x2 = (delta/2*(1+nodes))[None, :]
    alpha, kappa, cmu, cnuk = (0.45 if issuctionside else 0.30), .41, .09, 5.5
    ustar = uo*np.sqrt(cf/2)
    length = .085*delta*np.tanh(kappa*x2/(.085*delta))
    wake = 1-np.cos(pi*x2/delta)
    factor = uo/ustar-np.log(ustar*delta/kinvisc)/kappa-cnuk
    velocity = ustar*(np.log(ustar*x2/kinvisc)/kappa+cnuk+factor*.5*wake)
    gradient = ustar*(1/(kappa*x2)+factor*.5*pi/delta*np.sin(pi*x2/delta))
    ke = np.sqrt(pi)/length*.4213560764
    k1hat = k1/ke  # k3=0 in the upstream Pressure function
    nut = (length*kappa)**2*abs(gradient)
    kt = np.sqrt((nut*gradient)**2/cmu)
    ums = alpha*kt
    uc = .7*velocity
    ag = .05*uc/length
    phim = np.exp(-((omega-uc*k1)/ag)**2)/(ag*np.sqrt(pi))
    phi22 = 4/9/pi/ke**2*k1hat**2/(1+k1hat**2)**(7/3)
    inner = length*ums*gradient**2*phi22*phim*np.exp(-2*abs(k1)*x2)
    pressure = 4*airdens**2*(inner@weights)*delta/2
    return float(np.dot(weights, omega/spdsound/k1[:,0]*pressure)*(limits[1]-limits[0])/2)


def tblte_tno(u, theta, phi, d, r, cfall, d99all, edgevelall, p):
    """TBLTE_TNO. Includes upstream bandwidth, clipping, and Cf<=0 convention."""
    result = np.zeros((2,len(p.freqlist)))  # suction, pressure
    mach = u/p.spdsound
    directivity = directh_te(mach,theta,phi)
    ratio = 2**(1/3)
    for i, f in enumerate(p.freqlist):
        omega = twopi*f
        width = 2*omega*(np.sqrt(ratio)-1/np.sqrt(ratio))
        for side in (0,1):
            if cfall[side] > 0:
                answer = spl_integrate(omega,(0.,10*omega/u),side==0,mach,
                    p.spdsound,p.airdens,p.kinvisc,cfall,d99all,edgevelall)
                spectrum = d/(4*pi*r*r)*answer
                with np.errstate(divide='ignore'):
                    result[side,i] = 10*np.log10(spectrum*directivity/(2e-5)**2)+10*np.log10(width)
            result[side,i] = max(result[side,i], -100.)
    return result[1],result[0]


@dataclass
class BoundaryLayer:
    """Dimensional thicknesses (m); [suction, pressure] for each pair."""
    dstar: tuple
    d99: tuple
    cf: tuple = (0.,0.)
    edge_velocity_ratio: tuple = (1.,1.)


@dataclass
class BLTable:
    """AeroAcoustics ReadBLTables + BL_Param_Interp. Re is absolute, not millions."""
    aoa: np.ndarray
    reynolds: np.ndarray
    values: np.ndarray  # (nRe,nAoA,8): UeS,UeP,dstarS,dstarP,d99S,d99P,CfS,CfP

    @classmethod
    def read(cls, path):
        lines = Path(path).read_text(encoding='utf-8-sig').splitlines()
        nre, naoa = int(lines[2].split()[0]),int(lines[3].split()[0])
        tables, res, pos = [],[],4
        for _ in range(nre):
            while not lines[pos].strip() or lines[pos].lstrip().startswith('!'): pos+=1
            res.append(float(lines[pos].split()[0])*1e6)
            pos+=3
            table=np.array([[float(x) for x in line.split()[:9]] for line in lines[pos:pos+naoa]])
            if table.shape != (naoa,9): raise ValueError('Invalid boundary-layer table shape')
            tables.append(table); pos+=naoa
        tables=np.asarray(tables)
        aoa=(tables[0,:,0]+180)%360-180
        if np.any(np.diff(aoa)<=0) or np.any(np.diff(res)<=0):
            raise ValueError('BL table axes must be strictly increasing')
        if not np.allclose(tables[:,:,0],tables[0,:,0]) or not np.all(np.isfinite(tables)):
            raise ValueError('BL tables need equal AoA grids and finite values')
        return cls(aoa,np.asarray(res),tables[:,:,1:])

    def interpolate(self, alpha_deg, reynolds, chord):
        """Clamp isoparametric coordinates to [-1,1], as NWTC_Num does."""
        def bounds(grid,x):
            if len(grid)==1: return 0,0,0.
            low=int(np.clip(np.searchsorted(grid,x,side='right')-1,0,len(grid)-2))
            high=low+1
            return low,high,float(np.clip((x-grid[low])/(grid[high]-grid[low]),0.,1.))
        a,b,t=bounds(self.aoa,alpha_deg)
        c,d,s=bounds(self.reynolds,reynolds)
        v=(1-s)*((1-t)*self.values[c,a]+t*self.values[c,b])+s*((1-t)*self.values[d,a]+t*self.values[d,b])
        return BoundaryLayer(tuple(v[2:4]*chord),tuple(v[4:6]*chord),tuple(v[6:8]),tuple(v[:2]))


def observe(observer, aero_center, global_to_local, chord, airfoil_reference=(.25,0.)):
    """Calc_LE_Location_Array/CalcObserve for one node and one observer.

    Local x points toward suction, y along chord LE->TE, z root->tip.
    global_to_local is the 3x3 rotation supplied by the aerodynamic solver.
    Returns leading/trailing (distance, chordwise angle, spanwise angle), in m/deg.
    """
    rotation=np.asarray(global_to_local,dtype=float)
    centre=np.asarray(aero_center,dtype=float)
    if rotation.shape != (3,3) or not np.allclose(rotation@rotation.T,np.eye(3),atol=1e-8):
        raise ValueError('global_to_local must be an orthonormal 3x3 matrix')
    xr,yr=airfoil_reference
    output=[]
    for x in (0.,1.):
        offset=np.array([-yr,x-xr,0.])*chord
        location=offset@rotation+centre
        relative=rotation@(np.asarray(observer)-location)
        distance=max(aa_epsilon,float(np.linalg.norm(relative)))
        phi=np.arctan2(relative[0],relative[2])
        theta=np.arctan2(relative[2]*np.cos(phi)+relative[0]*np.sin(phi),relative[1])
        output.append((distance,float(np.degrees(theta)),float(np.degrees(phi))))
    return tuple(output)


def blade_elements(span, percentage=100.):
    """SetParameters node selection/lengths, including upstream edge conventions.

    Returns a zero-based start index and lengths. The source's startnode
    fallback can omit inner nodes even at 100%; retained for reproducibility.
    """
    span=np.asarray(span,dtype=float)
    if len(span)<1 or span[0]<0 or np.any(np.diff(span)<=0) or not 0<percentage<=100:
        raise ValueError('Increasing nonnegative span and 0<percentage<=100 required')
    n=len(span); start=max(1,n-1) # Fortran one-based index
    threshold=span[-1]*(1-percentage/100)
    for j in range(n-1,1,-1):
        if span[j-1] < threshold:
            start=j;break
    start=max(min(n,2),start)-1
    lengths=np.zeros(n)
    for j in range(start,n):
        if j==0: lengths[j]=span[j]
        elif j==n-1: lengths[j]=span[j]-span[j-1]
        else: lengths[j]=(span[j+1]-span[j-1])/2
    return start,lengths


class TurbulenceState:
    """AA_UpdateStates: supplied TI or 5 s regional moving buffers.

    Call after computing sound for the current step to reproduce the driver.
    Region sizes: 5 m radial, 60 degrees azimuth. Buffers are shared by blades.
    """
    def __init__(self,span,num_blades,dt,hub_height,method=1,ti=.1,avgv=8.):
        if dt<=0 or method not in (1,2): raise ValueError('Invalid dt/TICalcMeth')
        self.span=np.asarray(span);self.height=hub_height;self.method=method
        self.ti,self.avgv=ti,avgv
        self.n=max(int(np.floor(5/dt+.5)),1) # Fortran NINT
        self.radial=np.clip(np.ceil(self.span/5).astype(int),1,int(np.ceil(max(span)/5))+1)-1
        self.buffers=np.zeros((int(np.ceil(max(span)/5))+1,6,self.n))
        self.count=np.zeros(self.buffers.shape[:2],dtype=int)
        self.values=np.zeros((num_blades,len(span)))

    def update(self,vrel,inflow,leading_edge):
        vrel=np.asarray(vrel,dtype=float)
        if self.method==1:
            if np.any(vrel==0): raise ValueError('TI scaling requires nonzero Vrel')
            self.values[:]=self.ti*self.avgv/vrel
            return self.values.copy()
        for b in range(len(vrel)):
            for j in range(len(self.span)):
                z=leading_edge[b,j,2]-self.height
                if abs(z)<=max(abs(z),1.)*np.finfo(float).eps*50:
                    angular=0
                else:
                    angle=np.degrees(np.arctan2(leading_edge[b,j,1],z))%360
                    angular=int(np.clip(np.ceil(angle/60),1,6))-1
                radial=self.radial[j]
                self.count[radial,angular]+=1
                count=int(self.count[radial,angular]);buf=self.buffers[radial,angular]
                speed=float(np.linalg.norm(inflow[b,j]))
                if count<=self.n:
                    buf[count-1]=speed;self.values[b,j]=0.
                else:
                    buf[count%self.n]=speed  # preserve upstream MOD(counter,n)+1
                    mean=buf.mean()
                    self.values[b,j]=buf.std()/mean if abs(mean)>max(abs(mean),1.)*np.finfo(float).eps*50 else 0.
        return self.values.copy()


def section_spectrum(p, *, chord, speed, alpha_deg, span, distance=1.22,
                     theta=90.,phi=90.,leading_geometry=None,stall_deg=12.5,
                     boundary_layer=None,ti_section=None,te_thickness=.001,
                     te_angle=14.,thickness_1p=.02,thickness_10p=.12,is_tip=True):
    """CalcAeroAcousticsOutput for a section; shape (7,nFreq), in MECHANISMS order.

    Inactive sources are -inf (zero energy). Enabled routines keep the source's
    0 dB early-return and -100 dB floor conventions. No empirical formula fixes.
    """
    if not np.all(np.isfinite([chord,speed,alpha_deg,span,distance,theta,phi,stall_deg])):
        raise ValueError('Section inputs must be finite')
    if min(chord,span,distance)<=0 or speed<0 or speed>=p.spdsound:
        raise ValueError('Positive geometry and 0<=speed<sound speed required')
    speed=max(speed,.1)
    alpha=(alpha_deg+180)%360-180
    if p.tbltemod==2 and p.x_blmethod!=2:
        p=replace(p,x_blmethod=2)
    if p.x_blmethod==2 and boundary_layer is None:
        raise ValueError('BLMod=2 or TNO requires a boundary_layer or BL table')
    bl=boundary_layer or BoundaryLayer((0.,0.),(0.,0.))
    if p.x_blmethod==2 and (min(bl.dstar)<=0 or min(bl.d99)<=0):
        raise ValueError('Tabulated boundary-layer thicknesses must be positive')
    spl=np.full((7,len(p.freqlist)),-np.inf)
    common=(alpha,chord,speed,theta,phi,span,distance,p,bl.d99[1],bl.dstar[0],bl.dstar[1],stall_deg)
    if p.lammod==1 and p.itrip==0: spl[0]=lblvs(*common)
    if p.tbltemod:
        spl[1],spl[2],spl[3]=tblte(*common)
        if p.tbltemod==2:
            # CalcAeroAcousticsOutput explicitly overwrites EdgeVelVar with 1.
            spl[1],spl[2]=tblte_tno(speed,theta,phi,span,distance,bl.cf,bl.d99,(1.,1.),p)
    if p.bluntmod:
        if te_thickness<=0 or not 0<=te_angle<=14:
            raise ValueError('BPM bluntness requires h>0 and 0<=TE angle<=14 deg')
        spl[4]=blunt(alpha,chord,speed,theta,phi,span,distance,te_thickness,
                     te_angle,p,bl.d99[1],bl.dstar[0],bl.dstar[1],stall_deg)
    if p.tipmod and is_tip:
        spl[5]=tipnois(alpha,p.alprat,chord,speed,theta,phi,distance,p)
    if p.timod:
        intensity=p.ti*p.avgv/speed if ti_section is None else ti_section
        if intensity<0 or not np.isfinite(intensity): raise ValueError('Invalid sectional TI')
        r_le,t_le,p_le=leading_geometry or (distance,theta,phi)
        spl[6]=inflownoise(np.radians(alpha),chord,speed,t_le,p_le,span,r_le,intensity,p)
        if p.timod==2:
            spl[6]+=simple_guidati(speed,chord,thickness_10p,thickness_1p,p)+10.
    if p.aweighting: spl+=a_weighting(p.freqlist)
    return spl


def snapshot_spectrum(parameters, nodes, observers):
    """Pure-Python multi-node/multi-blade acoustic assembly.

    Each node supplies section arguments plus aero_center, global_to_local,
    and optional airfoil_reference. Return (observer,node,mechanism,frequency).
    Aerodynamic state is supplied, not inferred from a wind speed alone.
    """
    results=[]
    for observer in observers:
        spectra=[]
        for node in nodes:
            item=dict(node)
            centre=item.pop('aero_center');rotation=item.pop('global_to_local')
            reference=item.pop('airfoil_reference',(.25,0.))
            le,te=observe(observer,centre,rotation,item['chord'],reference)
            spectra.append(section_spectrum(parameters,**item,distance=te[0],theta=te[1],phi=te[2],leading_geometry=le))
        results.append(spectra)
    return np.asarray(results)


class AcousticDriver:
    """AeroAcoustics_Driver time loop + AA_CalcOutput/AA_UpdateStates.

    blades[b][j] is a node dictionary for snapshot_spectrum, plus inflow=[vx,vy,vz].
    span, is_tip and ti_section are supplied here using the upstream algorithms.
    Each step returns (observer,selected_node,mechanism,frequency), or None when
    AAStart/DT_AA suppresses output. State updates still occur after that step.
    """
    def __init__(self,p,span,num_blades,observers,dt=.1,aa_start=0.,
                 percentage=70.,hub_height=0.,ti_method=1):
        self.p=p;self.span=np.asarray(span,dtype=float);self.num_blades=num_blades
        self.observers=np.asarray(observers,dtype=float)
        self.dt=dt;self.start=aa_start;self.last_time=-np.inf
        self.first,self.lengths=blade_elements(self.span,percentage)
        self.state=TurbulenceState(span,num_blades,dt,hub_height,ti_method,p.ti,p.avgv)
        if self.observers.ndim!=2 or self.observers.shape[1]!=3 or len(self.observers)==0:
            raise ValueError('observers must have shape (n,3)')

    def step(self,time,blades):
        if not np.isfinite(time) or time<=self.last_time:
            raise ValueError('Time steps must be finite and strictly increasing')
        if len(blades)!=self.num_blades or any(len(b)!=len(self.span) for b in blades):
            raise ValueError('Blade/node counts differ from driver initialization')
        selected=[];shape=(self.num_blades,len(self.span))
        vrel=np.zeros(shape);inflow=np.zeros(shape+(3,));leading=np.zeros_like(inflow)
        for b,blade in enumerate(blades):
            for j,node in enumerate(blade):
                item=dict(node)
                inflow[b,j]=item.pop('inflow')
                vrel[b,j]=item['speed']
                reference=item.get('airfoil_reference',(.25,0.))
                local_le=np.array([-reference[1],-reference[0],0.])*item['chord']
                leading[b,j]=local_le@np.asarray(item['global_to_local'])+item['aero_center']
                if j>=self.first:
                    item.update(span=float(self.lengths[j]),is_tip=j==len(self.span)-1,
                                ti_section=float(self.state.values[b,j]))
                    selected.append(item)
        due=time>=self.start and (time+1e-10)%self.dt<1e-6
        output=snapshot_spectrum(self.p,selected,self.observers) if due else None
        self.state.update(vrel,inflow,leading)
        self.last_time=time
        return output


def read_observers(path):
    """Read AeroAcoustics observer file (count, header, then XYZ rows)."""
    lines=Path(path).read_text(encoding='utf-8-sig').splitlines()
    n=int(lines[0].split()[0])
    xyz=np.array([[float(x) for x in line.split()[:3]] for line in lines[2:2+n]])
    if xyz.shape!=(n,3) or n<1 or not np.all(np.isfinite(xyz)):
        raise ValueError('Invalid observer coordinates')
    return xyz


def read_aa_input(path):
    """Read named acoustic fields. Returns Parameters and scheduling/I/O metadata."""
    import re
    fields={}
    for line in Path(path).read_text(encoding='utf-8-sig').splitlines():
        match=re.match(r'''^\s*("[^"]*"|'[^']*'|\S+)\s+([A-Za-z]\w*)\b''',line)
        if match: fields[match[2].lower()]=match[1].strip('"\'')
    mapping={'spdsound':'spdsound','kinvisc':'kinvisc','airdens':'airdens',
        'lturb':'lturb','blmod':'x_blmethod','tripmod':'itrip','roundedtip':'round',
        'alprat':'alprat','timod':'timod','tbltemod':'tbltemod','lammod':'lammod',
        'tipmod':'tipmod','bluntmod':'bluntmod','aweighting':'aweighting','ti':'ti','avgv':'avgv'}
    ints={'x_blmethod','itrip','timod','tbltemod','lammod','tipmod','bluntmod'}
    values={}
    for key,dest in mapping.items():
        if key in fields:
            value=fields[key]
            if dest in ('round','aweighting'):
                if value.lower() not in ('true','false'): raise ValueError(f'Invalid {key}')
                values[dest]=value.lower()=='true'
            else: values[dest]=int(value) if dest in ints else float(value)
    return Parameters(**values), fields


def write_spectrum(output,p,spectrum,metadata):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    total=db_sum(spectrum,axis=0)
    with (output/'spectrum.csv').open('w',newline='',encoding='utf-8-sig') as stream:
        writer=csv.writer(stream);writer.writerow(('frequency_Hz',*MECHANISMS,'total_dB'))
        writer.writerows(zip(p.freqlist,*spectrum,total))
    finite_oaspl=lambda x: float(db_sum(x)) if np.any(np.isfinite(x)) else None
    summary={'source_commit':SOURCE_COMMIT,'weighting':'A' if p.aweighting else 'Z',
             'OASPL_dB':finite_oaspl(total),
             'mechanism_OASPL_dB':dict(zip(MECHANISMS,map(finite_oaspl,spectrum))),
             'case':metadata}
    (output/'summary.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8')
    return summary


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case',choices=['section','official-section'],default='section')
    parser.add_argument('--config',type=Path,help='JSON with parameters and section dictionaries')
    parser.add_argument('--tno',action='store_true',help='Demonstrate TNO instead of BPM TBL-TE')
    parser.add_argument('--a-weighting',action='store_true',help='Output A-weighted levels')
    parser.add_argument('--output',type=Path,default=Path('results'))
    parser.add_argument('--plot',action='store_true',help='Optional PNG; requires matplotlib')
    args=parser.parse_args()
    # Values explicitly present in upstream Aero_Tests; remaining inputs are stated assumptions.
    section=dict(chord=.22860,speed=63.920,alpha_deg=3.,span=.5090,distance=1.220,
                 theta=90.,phi=90.,stall_deg=12.5,te_thickness=.001,te_angle=14.)
    p=Parameters(itrip=0,lammod=1,tipmod=1,bluntmod=1)
    metadata={'name':'Aero_Tests section adaptation',
        'upstream_values':['chord','speed','alpha_deg','span','distance','theta','phi'],
        'assumptions':'Air properties, stall=12.5 deg, TI=.1, avgV=8, Lturb=40, h=.001 m, TE angle=14 deg; all mechanisms on; untripped BL.'}
    if args.case=='official-section':
        root=Path(__file__).resolve().parent/'reference/official_case'
        p,_=read_aa_input(root/'AeroAcousticsInput.dat')
        blade=np.loadtxt(root/'RotorSE_FAST_IEA_landBased_RWT_AeroDyn_blade.dat',skiprows=6)
        selected=int(np.flatnonzero(blade[:,6]==21)[0]) # AFID 21 maps to AF20
        chord=float(blade[selected,5]);length=float((blade[selected+1,0]-blade[selected-1,0])/2)
        table=BLTable.read(root/'Airfoils/AF20_BL.txt')
        coords=[]
        for line in (root/'Airfoils/AF20_Coords.txt').read_text().splitlines():
            try:
                parts=line.split()
                if len(parts)==2: coords.append([float(x) for x in parts])
            except ValueError: pass
        t1,t10=guidati_thickness(np.asarray(coords))
        section.update(chord=chord,span=length,speed=60.,distance=150.,alpha_deg=3.,stall_deg=8.500001,
            te_thickness=.01,te_angle=10.,
            boundary_layer=table.interpolate(3.,60.*chord/p.kinvisc,chord),thickness_1p=t1,thickness_10p=t10,is_tip=False)
        metadata={'name':'Official IEA AF20 section, prescribed aerodynamic state',
            'data':'Official AF20 coordinates/BL table and blade geometry',
            'assumptions':'Urel=60 m/s, AoA=3 deg, observer=150 m at theta=phi=90 deg; model switches follow official input.',
            'limitation':'An isolated official blade section; not a replay of the full OpenFAST rotor regression.'}
    if args.config:
        config=json.loads(args.config.read_text(encoding='utf-8-sig'))
        p=Parameters(**config['parameters']);section=config['section']
        if 'boundary_layer' in section: section['boundary_layer']=BoundaryLayer(**section['boundary_layer'])
        metadata={'name':str(args.config),'inputs':'User JSON'}
    if args.tno:
        p=replace(p,tbltemod=2,x_blmethod=2)
        if section.get('boundary_layer') is None:
            # Positive Cf/d99 values explicitly present in the Aero_Tests comments.
            # dstar for its extra separation term is generated by the BPM THICK routine.
            dp,ds_s,ds_p=thick(section['chord'],section['speed']*section['chord']/p.kinvisc,
                section['alpha_deg'],p,section['stall_deg'])
            section['boundary_layer']=BoundaryLayer((ds_s,ds_p),(.01105860,.007465830),(.000378576,.001984380))
        metadata['TNO_override']=True
    if args.a_weighting: p=replace(p,aweighting=True)
    metadata['parameters']=dict(asdict(p),freqlist=p.freqlist.tolist())
    metadata['section_inputs']=dict(section)
    if isinstance(section.get('boundary_layer'),BoundaryLayer):
        metadata['section_inputs']['boundary_layer']=asdict(section['boundary_layer'])
    spectrum=section_spectrum(p,**section)
    summary=write_spectrum(args.output,p,spectrum,metadata)
    if args.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig,ax=plt.subplots(figsize=(9,5.5),layout='constrained')
        for label,curve in zip(MECHANISMS,spectrum):
            if np.any(np.isfinite(curve)): ax.semilogx(p.freqlist,curve,label=label)
        ax.semilogx(p.freqlist,db_sum(spectrum,axis=0),'k',linewidth=2,label='total')
        ax.set(xlabel='1/3 octave centre frequency [Hz]',ylabel='SPL [dBA]' if p.aweighting else 'SPL [dB]',
            title=metadata['name'],ylim=(-20,max(80,float(np.max(spectrum))+5)))
        ax.grid(True,which='both',alpha=.2);ax.legend(ncol=4)
        fig.savefig(args.output/'spectrum.png',dpi=160);plt.close(fig)
    print(json.dumps(summary,ensure_ascii=False,indent=2))


def guidati_thickness(coords):
    """SetParameters nearest-coordinate sampling at 1%/10% chord.

    Input includes airfoil-reference row, then the TE->LE->TE contour.
    Initialize lower indices to the split node (upstream leaves them undefined
    if that node is already the nearest); no interpolation is introduced.
    """
    coords=np.asarray(coords,dtype=float)
    x,y=coords.T
    candidates=np.flatnonzero(np.diff(x[1:])>0)+2
    if len(candidates)==0: raise ValueError('Expected TE->LE->TE airfoil coordinates')
    split=int(candidates[0])
    return tuple(float(y[split+np.argmin(abs(x[split:]-target))]-y[np.argmin(abs(x[:split]-target))]) for target in (.01,.10))


if __name__=='__main__':
    main()
