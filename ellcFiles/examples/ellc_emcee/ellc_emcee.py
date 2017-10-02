#!/usr/bin/env python
 
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import numpy as np
from ellc import ldy,lc
import os
import emcee
from datetime import datetime
from astropy.table import Table, Column

def _lnlike(varpar_v, lcdata, varpar_l, fixpar_d, priors, ldy_, 
    n_spot_1, n_spot_2, grid_size, return_fit=False):

  try:
    r_1 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='r_1'])[0])]
  except:
    r_1 = fixpar_d['r_1']
  if (r_1 <= 0) or (r_1 > 1) : return -np.inf

  try:
    r_2 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='r_2'])[0])]
  except:
    r_2 = fixpar_d['r_2']
  if (r_2 <= 0) or (r_2 > 1) : return -np.inf

  try:
    sb2 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='sb2'])[0])]
  except:
    sb2 = fixpar_d['sb2']
  if sb2 < 0 : return -np.inf

  try:
    incl = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='incl'])[0])]
  except:
    incl = fixpar_d['incl']
  if (incl <= 0) or (incl > 90) : return -np.inf

  try:
    l_3 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='l_3'])[0])]
  except:
    l_3 = fixpar_d['l_3']

  try:
    t0 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='T_0'])[0])]
  except:
    t0 = fixpar_d['T_0']

  try:
    period = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='P'])[0])]
  except:
    period = fixpar_d['P']

  try:
    f_c = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='f_c'])[0])]
  except:
    f_c = fixpar_d['f_c']
  if (f_c <= -1) or (f_c >= 1): return -np.inf

  try:
    f_s = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='f_s'])[0])]
  except:
    f_s = fixpar_d['f_s']
  if (f_s <= -1) or (f_s >= 1): return -np.inf

  try:
    domdt = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='domdt'])[0])]
  except:
    domdt = fixpar_d['domdt']
  if (domdt <= -1) or (domdt >= 1): return -np.inf

  try:
    frot1 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='F_1'])[0])]
  except:
    frot1 = fixpar_d['F_1']

  try:
    frot2 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='F_2'])[0])]
  except:
    frot2 = fixpar_d['F_2']

  try:
    t1 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='Teff1'])[0])]
  except:
    t1 = fixpar_d['Teff1']

  try:
    g1 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='logg1'])[0])]
  except:
    g1 = fixpar_d['logg1']

  try:
    t2 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='Teff2'])[0])]
  except:
    t2 = fixpar_d['Teff2']

  try:
    g2 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='logg2'])[0])]
  except:
    g2 = fixpar_d['logg2']

  try:
    m = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='[M/H]'])[0])]
  except:
    m = fixpar_d['[M/H]']

  try:
    A_1 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='A_1'])[0])]
  except:
    A_1 = fixpar_d['A_1']
  if A_1 < 0 : return -np.inf

  try:
    A_2 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='A_2'])[0])]
  except:
    A_2 = fixpar_d['A_2']
  if A_2 < 0 : return -np.inf

  try:
    tilt = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='tilt'])[0])]
  except:
    tilt = fixpar_d['tilt']
  if tilt < 0 : return -np.inf

  spots_1 = None
  if n_spot_1 > 0:
    try:
      s_1_1 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='s_1_1'])[0])]
    except:
      s_1_1 = fixpar_d['s_1_1']
    if s_1_1 <0: return -np.inf
    try:
      l_1_1 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='l_1_1'])[0])]
    except:
      l_1_1 = fixpar_d['l_1_1']
    try:
      b_1_1 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='b_1_1'])[0])]
    except:
      b_1_1 = fixpar_d['b_1_1']
    try:
      f_1_1 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='f_1_1'])[0])]
    except:
      f_1_1 = fixpar_d['f_1_1']
    if f_1_1 <0: return -np.inf
    spots_1 = [[l_1_1],[b_1_1],[s_1_1],[f_1_1]]

  if n_spot_1 > 1:
    try:
      s_2_1 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='s_2_1'])[0])]
    except:
      s_2_1 = fixpar_d['s_2_1']
    if s_2_1 <0: return -np.inf
    try:
      l_2_1 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='l_2_1'])[0])]
    except:
      l_2_1 = fixpar_d['l_2_1']
    try:
      b_2_1 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='b_2_1'])[0])]
    except:
      b_2_1 = fixpar_d['b_2_1']
    try:
      f_2_1 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='f_2_1'])[0])]
    except:
      f_2_1 = fixpar_d['f_2_1']
    if f_2_1 <0: return -np.inf
    spots_1 = [[l_1_1,l_2_1],[b_1_1,b_2_1],[s_1_1,s_2_1],[f_1_1,f_2_1]]

  spots_2 = None
  if n_spot_2 > 0:
    try:
      s_1_2 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='s_1_2'])[0])]
    except:
      s_1_2 = fixpar_d['s_1_2']
    if s_1_2 <0: return -np.inf
    try:
      l_1_2 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='l_1_2'])[0])]
    except:
      l_1_2 = fixpar_d['l_1_2']
    try:
      b_1_2 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='b_1_2'])[0])]
    except:
      b_1_2 = fixpar_d['b_1_2']
    try:
      f_1_2 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='f_1_2'])[0])]
    except:
      f_1_2 = fixpar_d['f_1_2']
    if f_1_2 <0: return -np.inf
    spots_2 = [[l_1_2],[b_1_2],[s_1_2],[f_1_2]]

  if n_spot_2 > 1:
    try:
      s_2_2 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='s_2_2'])[0])]
    except:
      s_2_2 = fixpar_d['s_2_2']
    if s_2_2 <0: return -np.inf
    try:
      l_2_2 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='l_2_2'])[0])]
    except:
      l_2_2 = fixpar_d['l_2_2']
    try:
      b_2_2 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='b_2_2'])[0])]
    except:
      b_2_2 = fixpar_d['b_2_2']
    try:
      f_2_2 = varpar_v[(([i for i,s in enumerate(varpar_l) if s=='f_2_2'])[0])]
    except:
      f_2_2 = fixpar_d['f_2_2']
    if f_2_2 <0: return -np.inf
    spots_2 = [[l_1_2,l_2_2],[b_1_2,b_2_2],[s_1_2,s_2_2],[f_1_2,f_2_2]]



  a1,a2,a3,a4,y1 = ldy_(t1,g1,m)
  if np.isnan(y1):
    return -np.inf
  ldc_1 = [a1,a2,a3,a4]
  a1,a2,a3,a4,y2 = ldy_(t2,g2,m)
  if np.isnan(y2):
    return -np.inf
  ldc_2 = [a1,a2,a3,a4]

  heat_1 = A_1
  heat_2 = A_2

  q = fixpar_d['q']
  t_exp = fixpar_d['t_exp']
  if (t_exp <= 0) or (np.max(lcdata['flag']) < 2):
    t_exp = None

  if return_fit:
    n_int = np.array(lcdata['flag'])
    n_int[n_int<0] = 1 
    t = lcdata['time']
    f = lc(t, radius_1=r_1, radius_2=r_2, sbratio=sb2, incl=incl,
        light_3=l_3, t_zero=t0, period=period, q=q, f_c=f_c, f_s=f_s, 
        domdt=domdt,rotfac_1=frot1,rotfac_2=frot2,heat_1=heat_1,heat_2=heat_2,
        n_int=n_int, ld_1='claret', ldc_1=ldc_1, ld_2='claret', ldc_2=ldc_2,
        gdc_1=y1, gdc_2=y2, t_exp=t_exp, grid_1=grid_size, grid_2=grid_size,
        shape_1='roche', shape_2='roche',spots_1=spots_1,spots_2=spots_2)

    m = tilt*(t-np.min(t))-2.5*np.log10(f)
    m_mod = m[lcdata['flag'] >= 0]
    m_fit = (lcdata[lcdata['flag'] >= 0])['mag']
    e_fit = (lcdata[lcdata['flag'] >= 0])['e_mag']
    wt = 1.0/e_fit**2 
    zp = np.average( (m_fit-m_mod),weights=wt)
    return m+zp

  # Check for input value outside range set by priors
  for p in priors:

    if p['p_type'] == 'U':
      if p['p_var'] in varpar_l:
        v = varpar_v[(([i for i,s in enumerate(varpar_l) if s==p['p_var']])[0])]
      elif p['p_var'] == 'e':
        v = f_c**2 + f_s**2
      elif p['p_var'] == 'om':
        v = np.arctan2(f_s,f_c)
      elif p['p_var'] == 'rsum':
        v = r_1+r_2
      elif p['p_var'] == 'lrat':
        v = sb2*(r_2/r_1)**2
      elif p['p_var'] == 'k':
        v = r_2/r_1
      else:
        raise Exception('No such prior variable')

      if (v < p['p_par1']) or (v > p['p_par2']):
        return -np.inf

  t_fit = (lcdata[lcdata['flag'] >= 0])['time']
  n_int = (lcdata[lcdata['flag'] >= 0])['flag']
  f = lc(t_fit, radius_1=r_1, radius_2=r_2, sbratio=sb2, incl=incl, 
      light_3=l_3, t_zero=t0, period=period, q=q, f_c=f_c, f_s=f_s, domdt=domdt,
      rotfac_1=frot1, rotfac_2=frot2,heat_1=heat_1, heat_2=heat_2,
      n_int=n_int,ld_1='claret', ldc_1=ldc_1, ld_2='claret', ldc_2=ldc_2,
      gdc_1=y1, gdc_2=y2, t_exp=t_exp, grid_1=grid_size, grid_2=grid_size,
      shape_1='roche', shape_2='roche',spots_1=spots_1,spots_2=spots_2)


  if (True in np.isnan(f)) or np.min(f) <= 0 : return -np.inf

  mod_fit = tilt*(t_fit-np.min(t_fit))-2.5*np.log10(f)
  mag_fit = (lcdata[lcdata['flag'] >= 0])['mag']
  err_fit = (lcdata[lcdata['flag'] >= 0])['e_mag']
  wt = 1.0/err_fit**2
  zp = np.average( (mag_fit-mod_fit),weights=wt)
  mod_fit = mod_fit+zp
  lnlike = -0.5*np.sum(wt*(mod_fit-mag_fit)**2)

  for p in priors:
    if p['p_type'] == 'G':
      if p['p_var'] in varpar_l:
        v = varpar_v[(([i for i,s in enumerate(varpar_l) if s==p['p_var']])[0])]
      elif p['p_var'] == 'e':
        v = f_c**2 + f_s**2
      elif p['p_var'] == 'om':
        v = np.arctan2(f_s,f_c)
      elif p['p_var'] == 'rsum':
        v = r_1+r_2
      elif p['p_var'] == 'lrat':
        v = sb2*(r_2/r_1)**2
      elif p['p_var'] == 'k':
        v = r_2/r_1
      else:
        raise Exception('No such prior variable')

      lnlike = lnlike -0.5*((v-p['p_par1'])/(p['p_par2']))**2

  return lnlike
 
"""
 Light curve analysis using ellc binary star model and emcee 
"""

import argparse
from  astropy.io import fits
import numpy as np

def main():

  datetime_start = datetime.today()
  print("\n Start ellc_emcee.py at {:%c}\n".format(datetime_start))

  print(' Use ./ellc_emcee.py -h to see command line options.')
  print('\n Output file data format set by extension name, e.g. "chain.fits"')
  print(' See astropy.table docmentation for further details.')

  # Set up command line switched
  parser = argparse.ArgumentParser(
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    description='Light curve analysis using ellc binary star model and emcee'
  )
  parser.add_argument("--grid_size", "-g", 
    default='sparse',
    choices=['very_sparse', 'sparse', 'default', 'fine','very_fine'],
    help='Grid size'
  )
  parser.add_argument("--tolerant", "-t", 
    action="store_const",
    dest='tolerant',
    const=True,
    default=False,
    help='Tolerate invalid input'
  )
  parser.add_argument("--overwrite", "-f", 
    action="store_const",
    dest='overwrite',
    const=True,
    default=False,
    help='Force overwrite of existing output files.'
  )
  parser.add_argument("--log_file", "-l", 
    default='ellc_emcee.log',
    help='Log file name'
  )

  parser.add_argument("--chain_file", "-c", 
    default='chain.csv',
    help='Output file for chain data'
  )

  parser.add_argument("--model_file", "-m", 
    default='model.csv',
    help='Output file for the best fit model light curve.'
  )

  # Get command line options
  args = parser.parse_args()

  print(" Opening log file",args.log_file)
  if (args.overwrite):
    flags = os.O_CREAT | os.O_WRONLY
  else:
    flags = os.O_CREAT | os.O_EXCL | os.O_WRONLY

  try:
    file_handle = os.open(args.log_file,flags)
  except:
    raise

  f_log = os.fdopen(file_handle, 'w')

  print("\n Start ellc_emcee.py at {:%c}".format(datetime_start),file=f_log)

  print("\n Using ellc grid size",args.grid_size)
  print("\n Using ellc grid size",args.grid_size,file=f_log)

  print('\n Entering user options and parameter values.')
  print(' Blank lines and comments following "#" on input are ignored.')

  print("\n Light curve data file must be an ascii file with four columns of")
  print(" time, magnitude, error and flag, in that order, time in days. ")
  print(' Comment lines starting with "#" and additional columns are ignored.')
  print(' Flag is an integer used to set the following options.')
  print(' >= 0 to include observation in fit.')
  print(' > 1 to set no. of integration points to account for exposure time')
  print(' = 0 to interpolate the point from neighbouring values.')
  print(' < 0 to exclude point from fit, e.g., a point for plotting only.\n')

  prompt = "Enter light curve file name : "
  while True:
    line = raw_input(prompt)
    line = line.lstrip()
    if (len(line) == 0) and args.tolerant: continue   # Skip blank lines
    if line[0] == '#': continue   # Skip comment lines starting #
    line = (line.split('#'))[0]   # and ignore text after # on data lines

    lcfile = (line.split())[0]
    print("\n Attempting to read light curve file ",lcfile)
    print("\n Attempting to read light curve file ",lcfile,file=f_log)
    names = (str('time'),str('mag'),str('e_mag'),str('flag'))
    formats = ('f8','f4','f4','i4')
    try:
      lcdata = np.loadtxt(lcfile,dtype={'names': (names), 'formats': (formats)})
    except:
        if args.tolerant: 
          print ('Failed to read light curve data from file',lcfile)
          continue
        else:
          raise
    break
      
  print(' Succesfully read light curve file',lcfile)
  print(' Succesfully read light curve file',lcfile,file=f_log)
  n_data = len(lcdata)
  print(' Line  1:',lcdata[0])
  print(' Line {}:'.format(n_data),lcdata[n_data-1])
  print(' Line  1:',lcdata[0],file=f_log)
  print(' Line {}:'.format(n_data),lcdata[n_data-1],file=f_log)

  print('\n Bands: u,v,b,y (Stromgren)')
  print('        U, B, V, R, I, J, H, K (Johnson)')
  print('        S1, S2, S3, S4 (Spitzer IRAC)')
  print('        Kp (Kepler)')
  print('        C (CoRoT)')

  allowed_bands = ldy.list_bands()

  prompt = 'Enter two-character photometric band code for light curve :'
  while True:
    line = raw_input(prompt)
    line = line.lstrip()
    if (len(line) == 0) and args.tolerant: continue   # Skip blank lines
    if line[0] == '#': continue   # Skip comment lines starting #
    line = (line.split('#'))[0]   # and ignore text after # on data lines

    band = (line.split())[0]
    if band not in allowed_bands:
      if args.tolerant: 
        print ('Photometric band not available ',band)
        continue
      raise Exception ('Photometric band not available ',band)
    break

  print('\n Selected photometric band is',band)
  print('\n#  -- Start of user input block -- ',file=f_log)
  print(lcfile,'   #+ Light curve file ',file=f_log)
  print(band,'   #+ Photometric band ',file=f_log)

  def _get_param(prompt,default,lo_lim,hi_lim,integer=False):
    while True:
      if default is None:
        full_prompt = "{:s} : ".format(prompt)
      else:
        if integer:
          full_prompt = "{:s} [{:d}] : ".format(prompt,default)
        else:
          full_prompt = "{:s} [{:f}] : ".format(prompt,default)
      line = raw_input(full_prompt)
      if len(line) == 0 : 
        if default is None:
          continue
        else:
          print (default,'     #+ ',prompt,' (fixed, default)',file=f_log)
          return default,None
      if line[0] == '#': continue   # Skip comment lines starting #
      line = (line.split('#'))[0]   # and ignore text after # on data lines
      linedata = line.split()
      if len(linedata) == 0 : 
        if default is None:
          continue
        else:
          print (default,'     #+ ',prompt,' (fixed, default)',file=f_log)
          return default,None
      if integer:
        value = np.int(linedata[0])
      else:
        value = np.float(linedata[0])
      if (lo_lim is not None) and (value < lo_lim):
        if args.tolerant:
          print ('Input value exceeds lower limit',lo_lim)
          continue
        raise Exception('Input value exceeds lower limit',lo_lim)
      if (hi_lim is not None) and (value > hi_lim):
        if args.tolerant:
          print ('Input value exceeds upper limit',hi_lim)
          continue
        raise Exception('Input value exceeds upper limit',hi_lim)
      if len(linedata) == 1:
        print (value,'     #+ ',prompt,' (fixed)',file=f_log)
        return value,None
      else:
        error = np.float(linedata[1])
        print (value,error,'     #+ ',prompt,file=f_log)
        return value,error

  print('\n Now enter model parameters with their standard error estimate')
  print(' for free parameters, or value only for fixed parameters.')
  print(' T_eff, logg and [M/H] are only used to find limb darkening values.')
  print(' Hit enter to accept default (fixed) value where shown, e.g., [0].\n')

 # Could probably tidy this up by iterating over list of varaibles - oh well..

  varpar_v = []   # Values of free parameters
  varpar_e = []   # Errors on free parameters
  varpar_l = []   # List of free parameters
  fixpar_d = {}   # Dictionary of fixed parameters and values

  r_1,e_r1 = _get_param('r_1 = R_1/a',None,0.0001,1.0 )
  if e_r1 is None:
    fixpar_d['r_1'] = r_1
  else:
    varpar_v.append(r_1)
    varpar_e.append(e_r1)
    varpar_l.append('r_1')

  r_2,e_r2 = _get_param('r_2 = R_2/a',None,0.000001,100000)
  if e_r2 is None:
    fixpar_d['r_2'] = r_2
  else:
    varpar_v.append(r_2)
    varpar_e.append(e_r2)
    varpar_l.append('r_2')

  sbrat,e_sbrat = _get_param('s = S_2/S_1, surface brightness ratio',
      None,0.0,None)
  if e_sbrat is None:
    fixpar_d['sb2'] = sbrat
  else:
    varpar_v.append(sbrat)
    varpar_e.append(e_sbrat)
    varpar_l.append('sb2')

  incl,e_incl = _get_param('i = inclination (degrees)',None,0,90)
  if e_incl is None:
    fixpar_d['incl'] = incl
  else:
    varpar_v.append(incl)
    varpar_e.append(e_incl)
    varpar_l.append('incl')

  l_3,e_l_3 = _get_param('l_3 = third light',0,0,None)
  if e_l_3 is None:
    fixpar_d['l_3'] = l_3
  else:
    varpar_v.append(l_3)
    varpar_e.append(e_l_3)
    varpar_l.append('l_3')

  tzero,e_tzero = _get_param('T_0, time of mid-eclipse',None,None,None)
  if e_tzero is None:
    fixpar_d['T_0'] = tzero
  else:
    varpar_v.append(tzero)
    varpar_e.append(e_tzero)
    varpar_l.append('T_0')

  per,e_per = _get_param('P = orbital period (days)',None,0,None)
  if e_per is None:
    fixpar_d['P'] = per
  else:
    varpar_v.append(per)
    varpar_e.append(e_per)
    varpar_l.append('P')

  f_c,e_f_c = _get_param('f_c = sqrt(e).cos(omega)',0,-0.999,0.999)
  if e_f_c is None:
    fixpar_d['f_c'] = f_c
  else:
    varpar_v.append(f_c)
    varpar_e.append(e_f_c)
    varpar_l.append('f_c')

  f_s,e_f_s = _get_param('f_s = sqrt(e).sin(omega)',0,-0.999,0.999)
  if e_f_s is None:
    fixpar_d['f_s'] = f_s
  else:
    varpar_v.append(f_s)
    varpar_e.append(e_f_s)
    varpar_l.append('f_s')

  domdt,e_domdt = _get_param('domdt = apsidal motion rate [deg/cycle])',0,-1,1) 
  if e_domdt is None:
    fixpar_d['domdt'] = domdt
  else:
    varpar_v.append(domdt)
    varpar_e.append(e_domdt)
    varpar_l.append('domdt')

  frot1,e_frot1 = _get_param('F_1 = asynchronous rotation factor, star 1',
      1,None,None)
  if e_frot1 is None:
    fixpar_d['F_1'] = frot1
  else:
    varpar_v.append(frot1)
    varpar_e.append(e_frot1)
    varpar_l.append('F_1')

  frot2,e_frot2 = _get_param('F_2 = asynchronous rotation factor, star 2',
      1,None,None)
  if e_frot2 is None:
    fixpar_d['F_2'] = frot2
  else:
    varpar_v.append(frot2)
    varpar_e.append(e_frot2)
    varpar_l.append('F_2')

  teff1,e_teff1 = _get_param('Teff1 = T_eff,1',None,3500,50000)
  if e_teff1 is None:
    fixpar_d['Teff1'] = teff1
  else:
    varpar_v.append(teff1)
    varpar_e.append(e_teff1)
    varpar_l.append('Teff1')

  logg1,e_logg1 = _get_param('g1 = log(g_1)',4.0,0.0,5.0)
  if e_logg1 is None:
    fixpar_d['logg1'] = logg1
  else:
    varpar_v.append(logg1)
    varpar_e.append(e_logg1)
    varpar_l.append('logg1')

  teff2,e_teff2 = _get_param('Teff2 = T_eff,1',None,3500,50000)
  if e_teff2 is None:
    fixpar_d['Teff2'] = teff2
  else:
    varpar_v.append(teff2)
    varpar_e.append(e_teff2)
    varpar_l.append('Teff2')

  logg2,e_logg2 = _get_param('g2 = log(g_1)',4.0,0.0,5.0)
  if e_logg2 is None:
    fixpar_d['logg2'] = logg2
  else:
    varpar_v.append(logg2)
    varpar_e.append(e_logg2)
    varpar_l.append('logg2')

  M_H,e_M_H = _get_param('[M/H]',0.0,-5.0,1.0)
  if e_M_H is None:
    fixpar_d['[M/H]'] = M_H
  else:
    varpar_v.append(M_H)
    varpar_e.append(e_M_H)
    varpar_l.append('[M/H]')

  A_1,e_A_1 = _get_param('A_1 = reflection coefficient, star 1',0,0,None)
  if e_A_1 is None:
    fixpar_d['A_1'] = A_1
  else:
    varpar_v.append(A_1)
    varpar_e.append(e_A_1)
    varpar_l.append('A_1')

  A_2,e_A_2 = _get_param('A_2 = reflection coefficient, star 2',0,0,None)
  if e_A_2 is None:
    fixpar_d['A_2'] = A_2
  else:
    varpar_v.append(A_2)
    varpar_e.append(e_A_2)
    varpar_l.append('A_2')

  tilt,e_tilt = _get_param('tilt = linear trend d(mag)/dt',0,None,None)
  if e_tilt is None:
    fixpar_d['tilt'] = tilt
  else:
    varpar_v.append(tilt)
    varpar_e.append(e_tilt)
    varpar_l.append('tilt')

  print('\n Now enter fix parameters for the model. ')
  q,dummy = _get_param('Mass ratio q=m_2/m_1',1,0,None)  
  fixpar_d['q'] = q

  t_exp,dummy = _get_param('T_exp (seconds)',0,0,None)  
  fixpar_d['t_exp'] = t_exp/86400.

  print('\n Number of spots per star is limited to 2. ')

  allowed_prior_vars = ['e','om','r_1','r_2','rsum','k','sb2','incl','l_3',
      'lrat','T_0','P', 'f_c', 'f_s', 'domdt', 'F_1', 'F_2', 'A_1', 'A_2', 
      'Teff1', 'logg1', 'Teff2', 'logg2', '[M/H]'] 

  n_spot_1,dummy = _get_param('Number of spots on star 1',0,0,2,integer=True)  
  if n_spot_1 > 0:
    s_1_1,e_s_1_1 = _get_param('s_1_1 = size (degrees) of spot 1, star 1',
        None,0.001,60)  
    if e_s_1_1 is None:
      fixpar_d['s_1_1'] = s_1_1
    else:
      varpar_v.append(s_1_1)
      varpar_e.append(e_s_1_1)
      varpar_l.append('s_1_1')
      allowed_prior_vars.append('s_1_1')
    l_1_1,e_l_1_1 = _get_param('l_1_1 = longitude (degrees) of spot 1, star 1',
        None,-360,360)  
    if e_l_1_1 is None:
      fixpar_d['l_1_1'] = l_1_1
    else:
      varpar_v.append(l_1_1)
      varpar_e.append(e_l_1_1)
      varpar_l.append('l_1_1')
      allowed_prior_vars.append('l_1_1')
    b_1_1,e_b_1_1 = _get_param('b_1_1 = latitude (degrees) of spot 1, star 1',
        None,-90,90)
    if e_b_1_1 is None:
      fixpar_d['b_1_1'] = b_1_1
    else:
      varpar_v.append(b_1_1)
      varpar_e.append(e_b_1_1)
      varpar_l.append('b_1_1')
      allowed_prior_vars.append('b_1_1')
    f_1_1,e_f_1_1 = _get_param('f_1_1 = brightness factor of spot 1, star 1',
        None,0.0,None)  
    if e_f_1_1 is None:
      fixpar_d['f_1_1'] = f_1_1
    else:
      varpar_v.append(f_1_1)
      varpar_e.append(e_f_1_1)
      varpar_l.append('f_1_1')
      allowed_prior_vars.append('f_1_1')

  if n_spot_1 > 1:
    s_2_1,e_s_2_1 = _get_param('s_2_1 = size (degrees) of spot 1, star 1',
        None,0.001,60)  
    if e_s_2_1 is None:
      fixpar_d['s_2_1'] = s_2_1
    else:
      varpar_v.append(s_2_1)
      varpar_e.append(e_s_2_1)
      varpar_l.append('s_2_1')
      allowed_prior_vars.append('s_2_1')
    l_2_1,e_l_2_1 = _get_param('l_2_1 = longitude (degrees) of spot 1, star 1',
        None,-360,360)  
    if e_l_2_1 is None:
      fixpar_d['l_2_1'] = l_2_1
    else:
      varpar_v.append(l_2_1)
      varpar_e.append(e_l_2_1)
      varpar_l.append('l_2_1')
      allowed_prior_vars.append('l_2_1')
    b_2_1,e_b_2_1 = _get_param('b_2_1 = latitude (degrees) of spot 1, star 1',
        None,-90,90)
    if e_b_2_1 is None:
      fixpar_d['b_2_1'] = b_2_1
    else:
      varpar_v.append(b_2_1)
      varpar_e.append(e_b_2_1)
      varpar_l.append('b_2_1')
      allowed_prior_vars.append('b_2_1')
    f_2_1,e_f_2_1 = _get_param('f_2_1 = brightness factor of spot 1, star 1',
        None,0.0,None)  
    if e_f_2_1 is None:
      fixpar_d['f_2_1'] = f_2_1
    else:
      varpar_v.append(f_2_1)
      varpar_e.append(e_f_2_1)
      varpar_l.append('f_2_1')
      allowed_prior_vars.append('f_2_1')
   

  n_spot_2,dummy = _get_param('Number of spots on star 2',0,0,2,integer=True)  
  if n_spot_2 > 0:
    s_1_2,e_s_1_2 = _get_param('s_1_2 = size (degrees) of spot 1, star 2',
        None,0.001,60)  
    if e_s_1_2 is None:
      fixpar_d['s_1_2'] = s_1_2
    else:
      varpar_v.append(s_1_2)
      varpar_e.append(e_s_1_2)
      varpar_l.append('s_1_2')
      allowed_prior_vars.append('s_1_2')
    l_1_2,e_l_1_2 = _get_param('l_1_2 = longitude (degrees) of spot 1, star 2',
        None,-360,360)  
    if e_l_1_2 is None:
      fixpar_d['l_1_2'] = l_1_2
    else:
      varpar_v.append(l_1_2)
      varpar_e.append(e_l_1_2)
      varpar_l.append('l_1_2')
      allowed_prior_vars.append('l_1_2')
    b_1_2,e_b_1_2 = _get_param('b_1_2 = latitude (degrees) of spot 1, star 2',
        None,-90,90)
    if e_b_1_2 is None:
      fixpar_d['b_1_2'] = b_1_2
    else:
      varpar_v.append(b_1_2)
      varpar_e.append(e_b_1_2)
      varpar_l.append('b_1_2')
      allowed_prior_vars.append('b_1_2')
    f_1_2,e_f_1_2 = _get_param('f_1_2 = brightness factor of spot 1, star 2',
        None,0.0,None)  
    if e_f_1_2 is None:
      fixpar_d['f_1_2'] = f_1_2
    else:
      varpar_v.append(f_1_2)
      varpar_e.append(e_f_1_2)
      varpar_l.append('f_1_2')
      allowed_prior_vars.append('f_1_2')
   
  if n_spot_2 > 1:
    s_2_2,e_s_2_2 = _get_param('s_2_2 = size (degrees) of spot 1, star 2',
        None,0.001,60)  
    if e_s_2_2 is None:
      fixpar_d['s_2_2'] = s_2_2
    else:
      varpar_v.append(s_2_2)
      varpar_e.append(e_s_2_2)
      varpar_l.append('s_2_2')
      allowed_prior_vars.append('s_2_2')
    l_2_2,e_l_2_2 = _get_param('l_2_2 = longitude (degrees) of spot 1, star 2',
        None,-360,360)  
    if e_l_2_2 is None:
      fixpar_d['l_2_2'] = l_2_2
    else:
      varpar_v.append(l_2_2)
      varpar_e.append(e_l_2_2)
      varpar_l.append('l_2_2')
      allowed_prior_vars.append('l_2_2')
    b_2_2,e_b_2_2 = _get_param('b_2_2 = latitude (degrees) of spot 1, star 2',
        None,-90,90)
    if e_b_2_2 is None:
      fixpar_d['b_2_2'] = b_2_2
    else:
      varpar_v.append(b_2_2)
      varpar_e.append(e_b_2_2)
      varpar_l.append('b_2_2')
      allowed_prior_vars.append('b_2_2')
    f_2_2,e_f_2_2 = _get_param('f_2_2 = brightness factor of spot 1, star 2',
        None,0.0,None)  
    if e_f_2_2 is None:
      fixpar_d['f_2_2'] = f_2_2
    else:
      varpar_v.append(f_2_2)
      varpar_e.append(e_f_2_2)
      varpar_l.append('f_2_2')
      allowed_prior_vars.append('f_2_2')

  ndim = len(varpar_v)
  nwalkers,dummy = _get_param('Number of emcee walkers',
      ndim+1, 2*ndim, None, integer=True)  
  nsteps,dummy = _get_param('Number of steps in the emcee chain',
      1000, 1, None, integer=True)  
  nthreads,dummy = _get_param('Number of threads for emcee',
      1, 1, None, integer=True)  

  print('\n Now enter list of priors on parameters.\n')
  print(' In addition to parameters above (r_1, r_2, s, etc.), priors can also')
  print(' set on the following quantities:')
  print('  e = eccentricity')
  print('  om = longitude of periastron for star 1 (degrees)')
  print('  rsum = R_1/a + R_2/a')
  print('  k = R_2/R_1')
  print('  lrat = s*k^2 = luminosity ratio, L_2/L_1')

  print('\n Priors are entered using a parameter name followed by the type of ')
  print(' prior (U or G) and the parameters of the prior. ')
  print(' Prior types are as follows: ') 
  print('  U(lo_lim, hi_lim) - uniform between the limits lo_lim, hi_lim') 
  print('  G(mean, sigma) - Gaussian with specified mean and std. dev. = sigma')
  print(' Finish list of priors with "end". ')
  print(' e.g., ')
  print('  e G 0.15 0.05')
  print('  r_2 U 0.01 0.1') 
  print('  end \n')

  priors = np.zeros(0, dtype={
    'names':[str('p_var'), str('p_type'),str('p_par1'),str('p_par2')],
    'formats':['S5','S1','f8','f8']})

  allowed_prior_types = ['U','G']

  prompt = " Enter prior : var U/G lo_lim/mean hi_lim/sigma : "

  while True:
    line = raw_input(prompt)
    line = line.lstrip()
    line = (line.split('#'))[0]   # and ignore text after # on data lines
    if len(line) == 0 : continue  # Ignore blank lines/comment lines
    linedata = line.split()
    if (linedata[0]).lower() == "end": break
    if len(linedata) < 4:
      if args.tolerant: continue
      raise Exception ('Invalid entry ',line)
    p_var,p_type,p_str1,p_str2 = linedata[0:4]

    if p_var not in allowed_prior_vars:
      if args.tolerant:
        print ('No such variable',p_var)
        continue
      raise Exception ('No such variable',p_var)

    if p_type not in allowed_prior_types:
      if args.tolerant:
        print ('No such prior type',p_type)
        continue
      raise Exception ('No such prior type',p_type)

    try:
      p_par1 = np.float(p_str1)
    except:
      if args.tolerant:
        print('Failed to read floating point value ',p_str1)
        continue
      raise Exception('Failed to read floating point value ',p_str1)

    try:
      p_par2 = np.float(p_str2)
    except:
      if args.tolerant:
        print('Failed to read floating point value ',p_str2)
        continue
      raise Exception('Failed to read floating point value ',p_str2)

    
    if p_var in fixpar_d:
      print(' ### Warning - ignoring prior on fixed parameter : ',line)
      print(' ### Ignored prior on fixed parameter : ',line)
    elif p_var in priors['p_var']:
      print(' ### Warning - ignoring duplicate prior : ',line)
      print(' ### Ignored duplicate prior : ',line)
    else:
      priors = np.append(priors,
          np.array([(p_var,p_type,p_par1,p_par2)],dtype=priors.dtype))
      print('\n Including prior : ',line)
      print(line,'    #+  Prior',file=f_log)

  print('end    #+  End of priors',file=f_log)
  n_prior = len(priors)

  print('# --  End of user input block -- \n',file=f_log)
  print ('\n No. of free parameters = ',ndim)
  print ('\n No. of free parameters = ',ndim, file=f_log)
  print ('\n No. of priors = ',n_prior)
  print ('\n No. of priors = ',n_prior, file=f_log)
  t_fit = (lcdata[lcdata['flag'] >= 0])['time']
  n_fit = len(t_fit)
  n_df = n_fit + np.count_nonzero(priors['p_type'] == 'G') - ndim
  print(' Degrees of freedom (including Gaussian priors) =  ',n_df)
  print(' Degrees of freedom (including Gaussian priors) =  ',n_df,file=f_log)


# Limb- and gravity- darkening instance
  ldy_ = ldy.LimbGravityDarkeningCoeffs(band)

  # Initialize walkers so that none are out-of-bounds
  pos = []
  for i in range(nwalkers):
    lnlike_i = -np.inf
    while lnlike_i == -np.inf:
      pos_i = varpar_v + varpar_e*np.random.randn(ndim) 
      lnlike_i = _lnlike(pos_i,lcdata, varpar_l, fixpar_d,
          priors,ldy_,n_spot_1, n_spot_2,args.grid_size)
    pos.append(pos_i)

  sampler = emcee.EnsembleSampler(nwalkers, ndim, _lnlike,
      args=(lcdata, varpar_l, fixpar_d, priors,ldy_,n_spot_1, n_spot_2,
      args.grid_size),threads=nthreads)
  sampler.run_mcmc(pos, nsteps)
  af = sampler.acceptance_fraction
  print('\nMedian acceptance fraction =',np.median(af))
  print('\nMedian acceptance fraction =',np.median(af),file=f_log)

  t = Table(sampler.flatchain,names=varpar_l)
  t.add_column(Column(sampler.flatlnprobability,name='loglike'))
  indices = np.mgrid[0:nwalkers,0:nsteps]
  step = indices[1].flatten()
  walker = indices[0].flatten()
  t.add_column(Column(step,name='step'))
  t.add_column(Column(walker,name='walker'))
  t.write(args.chain_file)

  best_index = np.unravel_index(np.argmax(sampler.lnprobability),
      (nwalkers, nsteps))
  best_lnlike =  np.max(sampler.lnprobability)
  print('\n Best log(likelihood) = ',best_lnlike,' in walker ',best_index[0],
    ' at step ',best_index[1])

  best_pars = sampler.chain[best_index[0],best_index[1],:]

  lc_fit = _lnlike(best_pars, lcdata, varpar_l, fixpar_d, priors, ldy_,
      n_spot_1, n_spot_2, args.grid_size, return_fit=True)

  try:
    p = best_pars[(([i for i,s in enumerate(varpar_l) if s=='P'])[0])]
  except:
    p = fixpar_d['P']
  try:
    t0 = best_pars[(([i for i,s in enumerate(varpar_l) if s=='T_0'])[0])]
  except:
    t0 = fixpar_d['T_0']
  phase = ((lcdata['time'] - t0)/p) % 1

  t = Table(lcdata)
  t.add_column(Column(lc_fit,name='fit'))
  t.add_column(Column(phase,name='phase'))


  print("\n Writing best model fit to  ",args.model_file)
  print("\n Writing best model fit to  ",args.model_file,file=f_log)
  t.write(args.model_file)

  datetime_end = datetime.today()
  print("\n End ",__file__, "{:%c}\n".format(datetime_end))
  print("\n End ",__file__, "{:%c}\n".format(datetime_end),file=f_log)

#-------------------------------
if __name__ == "__main__":
    main()
