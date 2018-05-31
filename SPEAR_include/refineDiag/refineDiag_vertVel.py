#!/usr/bin/env python

import argparse
#import m6toolbox
import netCDF4 as nc
import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from verticalvelocity import calc_w_from_convergence
##-- RefineDiag Script for CMIP6
##
##   Variables we intend to provide in z-coordinates:
##
##     msftyyz    -> vmo  (ocean_z)  * both 0.25 and 0.5 resolutions
##     msftyzmpa  -> vhGM (ocean_z)  * applies only to 0.5 resolution
##     msftyzsmpa -> vhml (ocean_z)  * both 0.25 and 0.5 resolutions
##
##
##   Variables we intend to provide in rho-coordinates:
##   (potenital density referenced to 2000 m)
##
##     msftyrho    -> vmo
##     msftyrhompa -> vhGM           * applies only to 0.5 resolution
##
##
##   2-D variables we intent to provide:
##
##     hfy  ->  T_ady_2d + ndiff_tracer_trans_y_2d_T   * T_ady_2d needs to be converted to Watts (Cp = 3992.)
##                                                       ndiff_tracer_trans_y_2d_T already in Watts
##                                                       advective term in both 0.25 and 0.5 resolutions
##                                                       diffusive term only in 0.5 resolution
##
##     hfx  -> same recipie as above, expect for x-dimension
##     hfbasin -> summed line of hfy in each basin
##
##
##   Outstanding issues
##     1.) regirdding of vh, vhGM to rho-corrdinates
##     2.) vhGM and vhML units need to be in kg s-1
##     2.) save out RHO_0 and Cp somewhere in netCDF files to key off of
##
##
##   CMIP variables that will NOT be provided:
##
##     hfbasinpmadv, hfbasinpsmadv, hfbasinpmdiff, hfbasinpadv
##     (We advect the tracer with the residual mean velocity; individual terms cannot be diagnosed)
##
##     htovgyre, htovovrt, sltovgyre, sltovovrt
##     (Code for offline calculation not ready.)
##
##--

def run():
    parser = argparse.ArgumentParser(description='''CMIP6 RefineDiag Script for OM4''')
    parser.add_argument('infile', type=str, help='''Input file''')
    parser.add_argument('-b','--basinfile', type=str, default='', required=True, help='''File containing OM4 basin masks''')
    parser.add_argument('-o','--outfile', type=str, default=None, help='''Output file name''')
    parser.add_argument('-r','--refineDiagDir', type=str, default=None, help='''Path to refineDiagDir defined by FRE workflow)''')
    args = parser.parse_args()
    main(args)

def main(args):
    nc_misval = 1.e20
    print (args.basinfile)
    f_basin = nc.Dataset(args.basinfile)
    
    try:
      area_t = f_basin.variables['area_t'][:]
      wet    = f_basin.variables['wet'][:]
      do_areaT=True
      #print ("Found area in ",f_basin)
    except:
      do_areaT=False
      #print ("Could not find ",f_basin)

    #-- Read model data
    f_in = nc.Dataset(args.infile)

    #-- Read in existing dimensions from history netcdf file
    xh  = f_in.variables['xh']
    yh  = f_in.variables['yh']
    yq  = f_in.variables['yq']
    z_l = f_in.variables['z_l']
    z_i = f_in.variables['z_i']
    tax = f_in.variables['time']

    #-- Note: based on conversations with @adcroft, the overturning should be reported on the interfaces, z_i.
    #   Also, the nominal latitude is insufficient for the basin-average fields.  Based on the methods in
    #   meridional_overturning.py, the latitude dimension should be:
    #
    #   y = netCDF4.Dataset(cmdLineArgs.gridspec+'/ocean_hgrid.nc').variables['y'][::2,::2]
    #   yy = y[1:,:].max(axis=-1)+0*z
    #
    #   The quanity 'yy' above is numerically-equivalent to 'yq'

    #-- wmo
    if all(x in f_in.variables.keys() for x in ['umo', 'vmo']):
      varname = 'wmo'
      wmo = calc_w_from_convergence(f_in.variables['umo'], f_in.variables['vmo'])
      wmo[wmo.mask] = nc_misval
      wmo = np.ma.array(wmo,fill_value=nc_misval)
      wmo.long_name = 'Upward mass transport from resolved and parameterized advective transport'
      wmo.units = 'kg s-1'
      wmo.cell_methods = 'z_i:point xh:sum yh:sum time:mean'
      wmo.time_avg_info = 'average_T1,average_T2,average_DT'
      wmo.standard_name = 'upward_ocean_mass_transport'
      wmo.cell_measures = 'area:areacello'
      do_wmo = True
    else:
      do_wmo = False

    #-- wo
    if do_areaT:
      varname = 'wo'
      wo = calc_w_from_convergence(f_in.variables['umo'], f_in.variables['vmo'])
      # To avoid dividing nc_misval by area, set the masked out regions to 0.0
      wo[wo.mask] = 0.0
      wo = wo/(area_t*1035.0)
      #*wet
      wo[wo.mask] = nc_misval
      wo = np.ma.array(wo,fill_value=nc_misval)
      wo.long_name = 'Upward velocity from resolved and parameterized advective transport'
      wo.units = 'm s-1'
      wo.cell_methods = 'z_i:point xh:sum yh:sum time:mean'
      wo.time_avg_info = 'average_T1,average_T2,average_DT'
      wo.standard_name = 'upward_sea_water_velocity'
      wo.cell_measures = 'area:areacello'
      do_wo = True
    else:
      do_wo = False


    #-- Read time bounds
    nv = f_in.variables['nv']
    average_T1 = f_in.variables['average_T1']
    average_T2 = f_in.variables['average_T2']
    average_DT = f_in.variables['average_DT']
    time_bnds  = f_in.variables['time_bnds']

    if any([do_wmo, do_wo]):
      #-- Generate output filename
      if args.outfile is None:
        if hasattr(f_in,'filename'):
            args.outfile = f_in.filename
        else:
            args.outfile = os.path.basename(args.infile)
        args.outfile = args.outfile.split('.')
        args.outfile[-2] = args.outfile[-2]+'_refined'
        args.outfile = '.'.join(args.outfile)

      if args.refineDiagDir is not None:
        args.outfile = args.refineDiagDir + '/' + args.outfile

      #-- Write output file
      try:
          os.remove(args.outfile)
      except:
          pass

      if os.path.exists(args.outfile):
          raise IOError('Output netCDF file already exists.')
          exit(1)

      try:
        f_out     = nc.Dataset(args.outfile, 'w', format='NETCDF3_CLASSIC')
      except:
        f_out     = nc.Dataset(args.outfile, 'w')
      f_out.setncatts(f_in.__dict__)
      f_out.filename = os.path.basename(args.outfile)

      time_dim = f_out.createDimension('time', size=None)
      strlen_dim = f_out.createDimension('strlen', size=31)
      xh_dim  = f_out.createDimension('xh',  size=len(xh[:]))
      yh_dim  = f_out.createDimension('yh',  size=len(yh[:]))
      yq_dim  = f_out.createDimension('yq',  size=len(yq[:]))
      z_l_dim = f_out.createDimension('z_l', size=len(z_l[:]))
      z_i_dim = f_out.createDimension('z_i', size=len(z_i[:]))
      nv_dim  = f_out.createDimension('nv',  size=len(nv[:]))

      time_out = f_out.createVariable('time', np.float64, ('time'))
      xh_out   = f_out.createVariable('xh',   np.float64, ('xh'))
      yh_out   = f_out.createVariable('yh',   np.float64, ('yh'))
      yq_out   = f_out.createVariable('yq',   np.float64, ('yq'))
      z_l_out  = f_out.createVariable('z_l',  np.float64, ('z_l'))
      z_i_out  = f_out.createVariable('z_i',  np.float64, ('z_i'))
      nv_out  = f_out.createVariable('nv',  np.float64, ('nv'))

      if do_wmo:
        wmo_out = f_out.createVariable('wmo', np.float32, ('time', 'z_i', 'yh', 'xh'), fill_value=nc_misval)
        wmo_out.missing_value = nc_misval
        for k in wmo.__dict__.keys():
          if k[0] != '_': wmo_out.setncattr(k,wmo.__dict__[k])

      if do_wo:
        wo_out = f_out.createVariable('wo', np.float32, ('time', 'z_i', 'yh', 'xh'), fill_value=nc_misval)
        wo_out.missing_value = nc_misval
        for k in wo.__dict__.keys():
          if k[0] != '_': wo_out.setncattr(k,wo.__dict__[k])

      average_T1_out = f_out.createVariable('average_T1', np.float64, ('time'))
      average_T2_out = f_out.createVariable('average_T2', np.float64, ('time'))
      average_DT_out = f_out.createVariable('average_DT', np.float64, ('time'))
      time_bnds_out  = f_out.createVariable('time_bnds',  np.float64, ('time', 'nv'))

      time_out.setncatts(tax.__dict__)
      xh_out.setncatts(xh.__dict__)
      yh_out.setncatts(yh.__dict__)
      yq_out.setncatts(yq.__dict__)
      z_l_out.setncatts(z_l.__dict__)
      z_i_out.setncatts(z_i.__dict__)
      nv_out.setncatts(nv.__dict__)


      average_T1_out.setncatts(average_T1.__dict__)
      average_T2_out.setncatts(average_T2.__dict__)
      average_DT_out.setncatts(average_DT.__dict__)
      time_bnds_out.setncatts(time_bnds.__dict__)

      time_out[:] = np.array(tax[:])
      xh_out[:] = np.array(xh[:])
      yh_out[:] = np.array(yh[:])
      yq_out[:] = np.array(yq[:])
      z_l_out[:] = np.array(z_l[:])
      z_i_out[:] = np.array(z_i[:])
      nv_out[:] = np.array(nv[:])

      if do_wmo:        wmo_out[:] = np.ma.array(wmo[:])
      if do_wo:         wo_out[:]  = np.ma.array(wo[:])

      average_T1_out[:] = average_T1[:]
      average_T2_out[:] = average_T2[:]
      average_DT_out[:] = average_DT[:]
      time_bnds_out[:]  = time_bnds[:]


      f_out.close()
      exit(0)

    else:
      print('RefineDiag for ocean_month_rho2 yielded no output.')
      exit(1)



if __name__ == '__main__':
  run()
