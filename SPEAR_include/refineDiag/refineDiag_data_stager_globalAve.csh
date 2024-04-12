#------------------------------------------------------------------------------
#  refineDiag_data_stager_globalAve.csh
#
#  2014/05/07 JPK
#
#  DESCRIPTION:
#    This script serves two primary functions:
#
#    1.  It unpacks the history tar file to the /ptmp file system.  It allows
#        for more efficient post-processing when individual components are 
#        called by frepp.  (i.e. when the frepp "atmos_month" post-processing
#        script runs, frepp will copy only the unpacked "*atmos_month*" .nc 
#        files from /ptmp to the $work directory rather than the entire history
#        tar file.
#
#    2.  It performs a global annual average of all 3D variables (time, lat, lon)
#        and stores the values in a sqlite database that resides in a parallel
#        directory to the frepp scripts and stdout
#
#------------------------------------------------------------------------------
echo ""
echo ""
echo ""
echo "  ---------- begin refineDiag_data_stager.csh ----------  "
cd $work/$hsmdate
pwd
#-- This block is to run the refactored diagnostics and ingest them 
#-- into the mySQL database det_analysis currently running on Cobweb.
if ( $?tripleID  ) then
    if ( ${tripleID} != "" ) then
        python /home/fms/local/opt/fre-analysis/test/eem/code/detVitals/atmos_analysis.py -t ${tripleID} -y ${oname}
        python /home/fms/local/opt/fre-analysis/test/eem/code/detVitals/ocean_analysis.py -t ${tripleID} -y ${oname}
        python /home/fms/local/opt/fre-analysis/test/eem/code/detVitals/land_analysis.py -t ${tripleID} -y ${oname}
        python /home/fms/local/opt/fre-analysis/test/eem/code/detVitals/cobalt_analysis.py -t ${tripleID} -y ${oname}
        python /home/fms/local/opt/fre-analysis/test/eem/code/detVitals/detVitalsAverager.py -t ${tripleID} 
    else
        echo "tripleID for the experiment ${name} is not properly set... skipping detVitals diagnostics."
    endif
else
    echo "tripleID for the experiment ${name} is not properly set... skipping detVitals diagnostics."
endif
#-- Unload any previous versions of Python and load the system default
module unload python
module unload cdat
module load python
module load gcp
#-- Unpack gridSpec file.  Right now this hardcoded and this is bad practice.  
#   It would be much better to have the refineDiag script know about the gridSpec location
#   through an already populated FRE variable.  Will talk to Amy L. about alternatives.
#set gridSpecFile = "/archive/cjg/mdt/awg/input/grid/c96_GIS_025_grid_v20140327.tar"
#set gsArchRoot = `echo ${gridSpecFile} | rev | cut -f 2-100 -d '/' | rev`
#set gsBaseName = `basename ${gridSpecFile} | cut -f 1 -d '.'`
#hsmget -v -a ${gsArchRoot} -p /ptmp/$USER/${gsArchRoot} -w `pwd` ${gsBaseName}/\*
#-- Create a directory to house the sqlite database (if it does not already exist)
set localRoot = `echo $scriptName | rev | cut -f 4-100 -d '/' | rev`
if (! -d ${localRoot}/db) then 
  mkdir -p ${localRoot}/db
endif
set user = `whoami`
echo "oname is "${oname}
echo "hsmdate is "${hsmdate}
echo "gridspec is "${gridspec}
#-- Cat a Python script that performs the averages and writes to a copy of the DB
#   in case it is locked by another user

cat > global_atmos_ave.py <<EOF
import netCDF4 as nc
import numpy as np
import os
import pickle
import sqlite3
def getWebsiteVariablesDic():
    return pickle.load(open('/home/fms/local/opt/fre-analysis/test/eem/code/cm4_web_analysis/'+\
                            'etc/LM3_variable_dictionary.pkl', 'rb'))
def ncopen(file,action='exit'):
    if os.path.exists(file):
      return nc.Dataset(file)
    else:
      print('WARNING: Unable to open file '+file)
      if action == 'exit':
        exit(0)
      else:
        return None
def mask_latitude_bands(var,cellArea,geoLat,geoLon,region=None):
    if (region == 'tropics'):
      var = np.ma.masked_where(np.logical_or(geoLat < -30., geoLat > 30.),var)
      cellArea = np.ma.masked_where(np.logical_or(geoLat < -30., geoLat > 30.),cellArea)
    elif (region == 'nh'):
      var = np.ma.masked_where(np.less_equal(geoLat,30.),var)
      cellArea  = np.ma.masked_where(np.less_equal(geoLat,30.),cellArea)
    elif (region == 'sh'):
      var  = np.ma.masked_where(np.greater_equal(geoLat,-30.),var)
      cellArea  = np.ma.masked_where(np.greater_equal(geoLat,-30.),cellArea)
    elif (region == 'global'):
      var  = var
      cellArea = cellArea
    return var, cellArea
def area_mean(var,cellArea,geoLat,geoLon,cellFrac=None,soilFrac=None,region='global',varName=None,
              cellDepth=None, component=None):
    # Land-specific modifications
    if component == 'land':
        moduleDic = getWebsiteVariablesDic()
        # Read dictionary of keys
        if (varName in moduleDic.keys()):
          module = moduleDic[varName]
        elif (varName.lower() in moduleDic.keys()):
          module = moduleDic[varName.lower()]
        else:
          module = ''
        # Create a weighting factor
        if module == 'vegn':
          cellArea = cellArea*cellFrac*soilFrac
        else:
          cellArea = cellArea*cellFrac
        # Create a 3-D mask if needed
        if cellDepth is not None:
         if var.shape[0] == cellDepth.shape[0]:
           cellArea = np.tile(cellArea[None,:], (cellDepth.shape[0],1,1))
           geoLat = np.tile(geoLat[None,:], (cellDepth.shape[0],1,1))
           geoLon = np.tile(geoLon[None,:], (cellDepth.shape[0],1,1))
         else:
           print('Warning: inconsisent dimensions between varName and the cell depth axis.', \
                 var.shape[0], cellDepth.shape[0])
           null_result = np.ma.masked_where(True,0.)
           return null_result, null_result
        # Apply data mask to weighting mask
        cellArea.mask = var.mask
    var, cellArea = mask_latitude_bands(var,cellArea,geoLat,geoLon,region=region)
    #-- Land depth averaging and summation
    if cellDepth is not None:
      summed = np.ma.sum(var * cellArea * np.tile(cellDepth[:,None,None], (1,var.shape[1],var.shape[2])))
      var = np.ma.average(var,axis=0,weights=cellDepth)
      res = np.ma.sum(var*cellArea)/cellArea.sum()
      return res, summed
    else:
      res = np.ma.sum(var*cellArea)/cellArea.sum()
      return res, cellArea.sum()
def cube_sphere_aggregate(var,tiles):
    return np.ma.concatenate((tiles[0].variables[var][:], tiles[1].variables[var][:],\
                              tiles[2].variables[var][:], tiles[3].variables[var][:],\
                              tiles[4].variables[var][:], tiles[5].variables[var][:]),axis=-1)
def write_sqlite_data(sqlfile,varName,fYear,varmean=None,varsum=None,component=None):
    conn = sqlite3.connect(sqlfile)
    c = conn.cursor()
    if component == 'land':
      sql = 'create table if not exists '+varName+' (year integer primary key, sum float, avg float)'
    else:
      sql = 'create table if not exists '+varName+' (year integer primary key, value float)'
    sqlres = c.execute(sql)
    if component == 'land':
      sql = 'insert or replace into '+varName+' values('+fYear[:4]+','+str(varsum)+','+str(varmean)+')'
    else:
      sql = 'insert or replace into '+varName+' values('+fYear[:4]+','+str(varmean)+')'
    sqlres = c.execute(sql)
    conn.commit()
    c.close()
    conn.close()


def global_average_cubesphere(fYear, inputDir, outdir, label, history, ENSMEM):
    #import gmeantools
    import netCDF4 as nc
    import numpy as np
    import sqlite3
    import sys
    #fYearDir = inputDir + ENSMEM + "/history/" + fYear + ".nc/"+ fYear
    fYearDir = fYear
    gs_tiles = []
    #for tx in range(1,7): gs_tiles.append(gmeantools.ncopen(fYearDir + '.grid_spec.tile'+str(tx)+'.nc'))
    for tx in range(1,7): gs_tiles.append(ncopen(fYearDir + '.grid_spec.tile'+str(tx)+'.nc'))
    data_tiles = []
    #for tx in range(1,7): data_tiles.append(gmeantools.ncopen(fYearDir + '.'+history+'.tile'+str(tx)+'.nc'))
    #geoLat = gmeantools.cube_sphere_aggregate('grid_latt',gs_tiles)
    #geoLon = gmeantools.cube_sphere_aggregate('grid_lont',gs_tiles)
    #cellArea = gmeantools.cube_sphere_aggregate('area',gs_tiles)
    #for tx in range(1,7): data_tiles.append(ncopen(fYearDir + '.'+history+'.tile'+str(tx)+'.nc'))
    for tx in range(1,7): data_tiles.append(ncopen(fYearDir + '.'+history+'.tile'+str(tx)+'.nc'))
    geoLat = cube_sphere_aggregate('grid_latt',gs_tiles)
    geoLon = cube_sphere_aggregate('grid_lont',gs_tiles)
    cellArea = cube_sphere_aggregate('area',gs_tiles)
    for varName in data_tiles[0].variables.keys():
        if (len(data_tiles[0].variables[varName].shape) == 3):
            #var = gmeantools.cube_sphere_aggregate(varName,data_tiles)
            var = cube_sphere_aggregate(varName,data_tiles)
            var = np.ma.average(var,axis=0,weights=data_tiles[0].variables['average_DT'][:])
            for reg in ['global','tropics','nh','sh']:
                #result, _null = gmeantools.area_mean(var,cellArea,geoLat,geoLon,region=reg)
                result, _null = area_mean(var,cellArea,geoLat,geoLon,region=reg)
                #print(varName,reg,result)
                #gmeantools.
                write_sqlite_data(outdir+'/'+fYear+'.'+reg+'Ave'+label+'.db',varName,fYear[:4],result)


def global_average_land( fYear, inputDir, outdir, label, history, ENSMEM):
    #import gmeantools
    import numpy as np
    import netCDF4 as nc
    import pickle
    import re
    import sqlite3
    import sys
    #fYearDir = inputDir + ENSMEM + "/history/" + fYear + ".nc/"+ fYear
    fYearDir = fYear
    gs_tiles = []
    for tx in range(1,7): gs_tiles.append(ncopen(fYearDir + '.land_static.tile'+str(tx)+'.nc'))
    data_tiles = []
    for tx in range(1,7): data_tiles.append(ncopen(fYearDir + '.'+history+'.tile'+str(tx)+'.nc'))
    geoLat = cube_sphere_aggregate('geolat_t',gs_tiles)
    geoLon = cube_sphere_aggregate('geolon_t',gs_tiles)
    cellArea = cube_sphere_aggregate('land_area',gs_tiles)
    cellFrac = cube_sphere_aggregate('land_frac',gs_tiles)
    soilArea = cube_sphere_aggregate('soil_area',gs_tiles)
    soilFrac = np.ma.array(soilArea/(cellArea*cellFrac))
    depth = data_tiles[0].variables['zhalf_soil'][:]
    cellDepth = []
    for i in range(1,len(depth)):
        thickness = round((depth[i] - depth[i-1]),2)
        cellDepth.append(thickness)
    cellDepth = np.array(cellDepth)
    #for varName in data_tiles[0].variables.keys():
    for varName in ['LAI','frunf','hevap','levapv','npp','sens','soil_ice','evap_land','fsw','hprec','lwdn','precip','snow','soil_liq','flw','height','levapg','melt','runf','soil_T','transp']:

        varshape = data_tiles[0].variables[varName].shape
        if (len(varshape) >= 3):
            #var = gmeantools.cube_sphere_aggregate(varName,data_tiles)
            var = cube_sphere_aggregate(varName,data_tiles)
            var = np.ma.average(var,axis=0,weights=data_tiles[0].variables['average_DT'][:])
            if (len(varshape) == 3):
                for reg in ['global','tropics','nh','sh']:
                    #avg, wgt = gmeantools.area_mean(var,cellArea,geoLat,geoLon,cellFrac=cellFrac,soilFrac=soilFrac,\
                    avg, wgt = area_mean(var,cellArea,geoLat,geoLon,cellFrac=cellFrac,soilFrac=soilFrac,\
                                                   region=reg,varName=varName,component='land')
                    #print("3",varName,avg,avg*wgt)
                    if not hasattr(avg,'mask'):
                        #gmeantools.write_sqlite_data(outdir+'/'+fYear+'.'+reg+'Ave'+label+'.db',varName,fYear[:4],\
                        write_sqlite_data(outdir+'/'+fYear+'.'+reg+'Ave'+label+'.db',varName,fYear[:4],\
                                                     varmean=avg,varsum=avg*wgt,component='land')
        elif (len(varshape) == 4):
            if varshape[1] == cellDepth.shape[0]:
                for reg in ['global','tropics','nh','sh']:
                    #avg, summed = gmeantools.area_mean(var,cellArea,geoLat,geoLon,cellFrac=cellFrac,soilFrac=soilFrac,\
                    avg, summed = area_mean(var,cellArea,geoLat,geoLon,cellFrac=cellFrac,soilFrac=soilFrac,\
                                                       region=reg,varName=varName,cellDepth=cellDepth,component='land')
                    #print("4",varName,avg,summed)
                    #gmeantools.write_sqlite_data(outdir+'/'+fYear+'.'+reg+'Ave'+label+'.db',varName,fYear[:4],\
                    write_sqlite_data(outdir+'/'+fYear+'.'+reg+'Ave'+label+'.db',varName,fYear[:4],\
                                                 varmean=avg,varsum=summed,component='land')


def extract_ocean_scalar(fYear, inputDir, outdir, ENSMEM):
    #import gmeantools
    import netCDF4 as nc
    import numpy as np
    import sqlite3
    import sys
    #fYearDir = inputDir + ENSMEM + "/history/" + fYear + ".nc/"+ fYear
    fYearDir = fYear
    fdata = ncopen(fYearDir + '.ocean_scalar_month.nc')
    ignoreList = ['time_bounds', 'time_bnds', 'average_T2', 'average_T1', 'average_DT']
    varDict = fdata.variables.keys()
    varDict = list(set(varDict) - set(ignoreList))
    for varName in varDict:
        if len(fdata.variables[varName].shape) == 2:
            result = np.ma.average(fdata.variables[varName][:,:],axis=0,weights=fdata.variables['average_DT'][:]).sum()
            write_sqlite_data(outdir+'/'+fYear+'.globalAveOcean.db',varName,fYear[:4],result)


global_average_cubesphere(str('${oname}'), '.', '.', "Atmos", "atmos_month", '')
global_average_cubesphere(str('${oname}'), '.', '.', "AtmosAer", "atmos_month_aer", '')

global_average_land( str('${oname}'), '.', ".", "Land", "land_month", '')

#Move to global_oean_ave.py
#extract_ocean_scalar( str('${oname}'), '.','.', '')

EOF

cat > global_ocean_ave.py <<EOF
import netCDF4 as nc
import numpy as np
import os
import pickle
import sqlite3

def getWebsiteVariablesDic():
    return pickle.load(open('/home/fms/local/opt/fre-analysis/test/eem/code/cm4_web_analysis/'+\
                            'etc/LM3_variable_dictionary.pkl', 'rb'))

def ncopen(file,action='exit'):
    if os.path.exists(file):
      return nc.Dataset(file)
    else:
      print('WARNING: Unable to open file '+file)
      if action == 'exit':
        exit(0)
      else:
        return None

def mask_latitude_bands(var,cellArea,geoLat,geoLon,region=None):
    if (region == 'tropics'):
      var = np.ma.masked_where(np.logical_or(geoLat < -30., geoLat > 30.),var)
      cellArea = np.ma.masked_where(np.logical_or(geoLat < -30., geoLat > 30.),cellArea)
    elif (region == 'nh'):
      var = np.ma.masked_where(np.less_equal(geoLat,30.),var)
      cellArea  = np.ma.masked_where(np.less_equal(geoLat,30.),cellArea)
    elif (region == 'sh'):
      var  = np.ma.masked_where(np.greater_equal(geoLat,-30.),var)
      cellArea  = np.ma.masked_where(np.greater_equal(geoLat,-30.),cellArea)
    elif (region == 'global'):
      var  = var
      cellArea = cellArea
    return var, cellArea

def area_mean(var,cellArea,geoLat,geoLon,cellFrac=None,soilFrac=None,region='global',varName=None,
              cellDepth=None, component=None):
    # Land-specific modifications
    if component == 'land':
        moduleDic = getWebsiteVariablesDic()
        # Read dictionary of keys
        if (varName in moduleDic.keys()):
          module = moduleDic[varName]
        elif (varName.lower() in moduleDic.keys()):
          module = moduleDic[varName.lower()]
        else:
          module = ''
        # Create a weighting factor
        if module == 'vegn':
          cellArea = cellArea*cellFrac*soilFrac
        else:
          cellArea = cellArea*cellFrac
        # Create a 3-D mask if needed
        if cellDepth is not None:
         if var.shape[0] == cellDepth.shape[0]:
           cellArea = np.tile(cellArea[None,:], (cellDepth.shape[0],1,1))
           geoLat = np.tile(geoLat[None,:], (cellDepth.shape[0],1,1))
           geoLon = np.tile(geoLon[None,:], (cellDepth.shape[0],1,1))
         else:
           print('Warning: inconsisent dimensions between varName and the cell depth axis.', \
                 var.shape[0], cellDepth.shape[0])
           null_result = np.ma.masked_where(True,0.)
           return null_result, null_result
        # Apply data mask to weighting mask
        cellArea.mask = var.mask
    var, cellArea = mask_latitude_bands(var,cellArea,geoLat,geoLon,region=region)
    #-- Land depth averaging and summation
    if cellDepth is not None:
      summed = np.ma.sum(var * cellArea * np.tile(cellDepth[:,None,None], (1,var.shape[1],var.shape[2])))
      var = np.ma.average(var,axis=0,weights=cellDepth)
      res = np.ma.sum(var*cellArea)/cellArea.sum()
      return res, summed
    else:
      res = np.ma.sum(var*cellArea)/cellArea.sum()
      return res, cellArea.sum()

def cube_sphere_aggregate(var,tiles):
    return np.ma.concatenate((tiles[0].variables[var][:], tiles[1].variables[var][:],\
                              tiles[2].variables[var][:], tiles[3].variables[var][:],\
                              tiles[4].variables[var][:], tiles[5].variables[var][:]),axis=-1)


def write_sqlite_data(sqlfile,varName,fYear,varmean=None,varsum=None,component=None):
    conn = sqlite3.connect(sqlfile)
    c = conn.cursor()
    if component == 'land':
      sql = 'create table if not exists '+varName+' (year integer primary key, sum float, avg float)'
    else:
      sql = 'create table if not exists '+varName+' (year integer primary key, value float)'
    sqlres = c.execute(sql)
    if component == 'land':
      sql = 'insert or replace into '+varName+' values('+fYear[:4]+','+str(varsum)+','+str(varmean)+')'
    else:
      sql = 'insert or replace into '+varName+' values('+fYear[:4]+','+str(varmean)+')'
    sqlres = c.execute(sql)
    conn.commit()
    c.close()
    conn.close()


import netCDF4 as nc
import tarfile
import numpy
import sqlite3
from scipy.io import netcdf
#import gmeantools

def ice9it(i, j, depth, minD=0.):
    wetMask = 0*depth      
    (nj,ni) = wetMask.shape
    stack = set()
    stack.add( (j,i) )
    while stack:
        (j,i) = stack.pop()
        if wetMask[j,i] or depth[j,i] <= minD: continue
        wetMask[j,i] = 1
  
        if i>0: stack.add( (j,i-1) )
        else: stack.add( (j,ni-1) ) # Periodic beyond i=0
        if i<ni-1: stack.add( (j,i+1) )
        else: stack.add( (j,0) ) # Periodic beyond i=ni-1
        
        if j>0: stack.add((j-1,i))
        
        if j<nj-1: stack.add( (j+1,i) )
        else: stack.add( (j,ni-1-i) ) # Tri-polar fold beyond j=nj-1
    return wetMask

def ice9(x, y, depth, xy0):
    #ji = nearestJI(x, y, xy0)
    ji = nearestJI(x, y, xy0[0],xy0[1])
    return ice9it(ji[1], ji[0], depth)

def nearestJI(x, y, x0, y0):
    return numpy.unravel_index( ((x-x0)**2 + (y-y0)**2).argmin() , x.shape)

def southOf(x, y, xy0, xy1):
    x0 = xy0[0]; y0 = xy0[1]; x1 = xy1[0]; y1 = xy1[1]
    dx = x1 - x0; dy = y1 - y0
    Y = (x-x0)*dy - (y-y0)*dx
    Y[Y>=0] = 1; Y[Y<=0] = 0
    return Y


def MOCpsi(vh, vmsk=None, yq=None):
    shapeMask = list(vmsk.shape);
    shape = list(vh.shape); 
    #print(shape, shapeMask, type(vmsk))
    shape[-3] += 1
    #print(shape, shapeMask, type(vmsk))
    if len(shape)==3:
        #print("len shape = 3")
        if vmsk is not None:
            if shape[1]==shapeMask[0]+1:
                # Symmetric grid (ignore southern most latitude)
                #print "vh_=vh[:,1:,:] shape1"
                vh_=vh[:,1:,:]
                yq_=yq[1:]
            elif shape[2]==shapeMask[1]+1:  
                vh_=vh[:,:,1:]
                yq_=yq[1:]
                #print "vh_=vh[:,:,1:] shape2"
            else:
                vh_=vh
                yq_=yq
        shapevh = list(vh_.shape); 
        #print( "shape vh_",shapevh)
        #shapevh[-3] +=1
        psi = numpy.zeros(shapevh[:-1])
        for k in range(shapevh[-3]-1,0,-1):
            #print 'k=',k
            if vmsk is not None: 
                psi[k-1,:] = psi[k,:] - vh_[k-1].sum(axis=-1)
            else: 
                psi[k-1,:] = psi[k,:] - (vmsk*vh_[k-1]).sum(axis=-1)
    else:
        if vmsk is not None:
            if shape[2]==shapeMask[0]+1:  
                vh_=vh[:,:,1:,:]
                yq_=yq[1:]
            elif shape[3]==shapeMask[1]+1:  
                vh_=vh[:,:,:,1:]
                yq_=yq[1:]
            else:
                vh_=vh
        for n in range(shape[0]):
            for k in range(shape[-3]-1,0,-1):
                if vmsk==None: psi[n,k-1,:] = psi[n,k,:] - vh_[n,k-1].sum(axis=-1)
                else: psi[n,k-1,:] = psi[n,k,:] - (vmsk*vh_[n,k-1]).sum(axis=-1)
    return psi,yq_




def grid_info_from_gridspec(gridSpec):
    #-- Get grid info from gridspec file
    gsFile = gridSpec #"${gridspec}"
    TF = tarfile.open(gsFile,'r')
    member = [m for m in TF.getmembers() if 'ocean_hgrid' in m.name][0]
    nc = netcdf.netcdf_file(TF.extractfile(member),'r')
    x = nc.variables['x'][1::2,1::2]
    y = nc.variables['y'][1::2,1::2]
    member = [m for m in TF.getmembers() if 'topog' in m.name][0]
    nc = netcdf.netcdf_file(TF.extractfile(member),'r')
    depth = nc.variables['depth'][:]
    #print 'Generating global wet mask ...',
    wet = ice9(x, y, depth, (0,-35)) # All ocean points seeded from South Atlantic
    #print 'done.'
    code = 0*wet
    #print ('Finding Cape of Good Hope ...')
    tmp = 1 - wet; tmp[x<-30] = 0
    tmp = ice9(x, y, tmp, (20,-30.))
    yCGH = (tmp*y).min()
    #print( 'done.', yCGH)
    #print ('Finding Melbourne ...')
    tmp = 1 - wet; tmp[x>-180] = 0
    tmp = ice9(x, y, tmp, (-220,-25.))
    yMel = (tmp*y).min()
    #print ('done.', yMel)
    #print ('Processing Persian Gulf ...')
    tmp = wet*( 1-southOf(x, y, (55.,23.), (56.5,27.)) )
    tmp = ice9(x, y, tmp, (53.,25.))
    code[tmp>0] = 11
    wet = wet - tmp # Removed named points
    #print ('Processing Red Sea ...')
    tmp = wet*( 1-southOf(x, y, (40.,11.), (45.,13.)) )
    tmp = ice9(x, y, tmp, (40.,18.))
    code[tmp>0] = 10
    wet = wet - tmp # Removed named points
    #print 'Processing Black Sea ...'
    tmp = wet*( 1-southOf(x, y, (26.,42.), (32.,40.)) )
    tmp = ice9(x, y, tmp, (32.,43.))
    code[tmp>0] = 7
    wet = wet - tmp # Removed named points
    #print 'Processing Mediterranean ...'
    tmp = wet*( southOf(x, y, (-5.7,35.5), (-5.7,36.5)) )
    tmp = ice9(x, y, tmp, (4.,38.))
    code[tmp>0] = 6
    wet = wet - tmp # Removed named points
    #print 'Processing Baltic ...'
    tmp = wet*( southOf(x, y, (8.6,56.), (8.6,60.)) )
    tmp = ice9(x, y, tmp, (10.,58.))
    code[tmp>0] = 9
    wet = wet - tmp # Removed named points
    #print 'Processing Hudson Bay ...'
    tmp = wet*( 
               ( 1-(1-southOf(x, y, (-95.,66.), (-83.5,67.5)))
                  *(1-southOf(x, y, (-83.5,67.5), (-84.,71.))) 
               )*( 1-southOf(x, y, (-70.,58.), (-70.,65.)) ) )
    tmp = ice9(x, y, tmp, (-85.,60.))
    code[tmp>0] = 8
    wet = wet - tmp # Removed named points
    #print 'Processing Arctic ...'
    tmp = wet*( 
              (1-southOf(x, y, (-171.,66.), (-166.,65.5))) * (1-southOf(x, y, (-64.,66.4), (-50.,68.5))) # Lab Sea
         +    southOf(x, y, (-50.,0.), (-50.,90.)) * (1- southOf(x, y, (0.,65.5), (360.,65.5))  ) # Denmark Strait
         +    southOf(x, y, (-18.,0.), (-18.,65.)) * (1- southOf(x, y, (0.,64.9), (360.,64.9))  ) # Iceland-Sweden
         +    southOf(x, y, (20.,0.), (20.,90.)) # Barents Sea
         +    (1-southOf(x, y, (-280.,55.), (-200.,65.)))
              )
    tmp = ice9(x, y, tmp, (0.,85.))
    code[tmp>0] = 4
    wet = wet - tmp # Removed named points
    #print 'Processing Pacific ...'
    tmp = wet*( (1-southOf(x, y, (0.,yMel), (360.,yMel)))
               -southOf(x, y, (-257,1), (-257,0))*southOf(x, y, (0,3), (1,3))
               -southOf(x, y, (-254.25,1), (-254.25,0))*southOf(x, y, (0,-5), (1,-5)) 
               -southOf(x, y, (-243.7,1), (-243.7,0))*southOf(x, y, (0,-8.4), (1,-8.4)) 
               -southOf(x, y, (-234.5,1), (-234.5,0))*southOf(x, y, (0,-8.9), (1,-8.9)) 
              )
    tmp = ice9(x, y, tmp, (-150.,0.))
    code[tmp>0] = 3
    wet = wet - tmp # Removed named points
    #print 'Processing Atlantic ...'
    tmp = wet*(1-southOf(x, y, (0.,yCGH), (360.,yCGH)))
    tmp = ice9(x, y, tmp, (-20.,0.))
    code[tmp>0] = 2
    wet = wet - tmp # Removed named points
    #print 'Processing Indian ...'
    tmp = wet*(1-southOf(x, y, (0.,yCGH), (360.,yCGH)))
    tmp = ice9(x, y, tmp, (55.,0.))
    code[tmp>0] = 5
    wet = wet - tmp # Removed named points
    #print 'Processing Southern Ocean ...'
    tmp = ice9(x, y, wet, (0.,-55.))
    code[tmp>0] = 1
    wet = wet - tmp # Removed named points
    code[wet>0] = -9
    (j,i) = numpy.unravel_index( wet.argmax(), x.shape)
    if j:
      print ('There are leftover points unassigned to a basin code')
      while j:
        print (x[j,i],y[j,i],[j,i])
        wet[j,i]=0
        (j,i) = numpy.unravel_index( wet.argmax(), x.shape)
    else: print ('All points assigned a basin code')
    #-- Define atlantic/arctic mask
    atlmask = numpy.where(numpy.logical_or(code==2,code==4),1.,0.)
    return atlmask
    


def Read_VMO_VH(fYearDir):
    import xarray as xr
    #-- Attempt to read vmo and/or vh
    doVmo = False
    doVh = False
    vmo = None
    zl = None
    yq_vmo = None
    yq_vh = None
    try:
        vhFile = xr.open_dataset(fYearDir + '.ocean_z_month.nc')
        if 'vmo' in vhFile.variables.keys():
            vmo=vhFile['vmo'].resample(time="1Y").mean()*1e-9
            #print("annual shape is ",annual.shape)
            #vmo  = annual.filled(0)*1e-9
            #print("vmo shape is ",vmo.shape)
            zl  = vhFile.variables['z_l'][:]
            yq_vmo  = vhFile.variables['yq'][:]
            doVmo = True
    except: doVh = False
    return vmo[0],zl,yq_vmo,yq_vh,doVmo,doVh
    

def extract_ocean_scalar(fYear, inputDir, outdir, ENSMEM):
    #import gmeantools
    import netCDF4 as nc
    import numpy as np
    import sqlite3
    import sys
    #fYearDir = inputDir + ENSMEM + "/history/" + fYear + ".nc/"+ fYear
    fYearDir = fYear
    fdata = ncopen(fYearDir + '.ocean_scalar_month.nc')
    ignoreList = ['time_bounds', 'time_bnds', 'average_T2', 'average_T1', 'average_DT']
    varDict = fdata.variables.keys()
    varDict = list(set(varDict) - set(ignoreList))
    for varName in varDict:
        if len(fdata.variables[varName].shape) == 2:
            result = np.ma.average(fdata.variables[varName][:,:],axis=0,weights=fdata.variables['average_DT'][:]).sum()
            #print(outdir+'/'+fYear+'.globalAveOcean.db',varName,fYear[:4],result)
            write_sqlite_data(outdir+'/'+fYear+'.globalAveOcean.db',varName,fYear[:4],result)
    #Now calculate AMOC 
    vmo,zl,yq_vmo,yq_vh,doVmo,doVh = Read_VMO_VH(fYearDir)
    amocDict = {}
    if doVmo == True:
        psi,yq_vmo_ = MOCpsi(vmo,vmsk=atlmask,yq=yq_vmo)
        #print  ("shape of zl",list(zl.shape))
        #print  ("shape of MOCpsi",list(psi.shape))
        #print  ("shape of yq_vmo",list(yq_vmo.shape))
        #print  ("shape of yq_vmo_",list(yq_vmo_.shape))
        maxsfn = numpy.max(psi[numpy.logical_and(zl>500,zl<2500)][:,numpy.greater_equal(yq_vmo_,20)])
        #print('AMOC vmo = %s' % maxsfn)
        amocDict['amoc_vmo'] = maxsfn
    if doVh == True:
        psi,yq_vh_ = MOCpsi(vh,vmsk=atlmask,yq=yq_vh)
        maxsfn = numpy.max(psi[numpy.logical_and(zt>500,zt<2500)][:,numpy.greater_equal(yq_vh_,20)])
        #print('AMOC vh = %s' % maxsfn)
        amocDict['amoc_vh'] = maxsfn
    for varName in amocDict:
        result= amocDict[varName]
        #print(outdir+'/'+fYear+'.globalAveOcean.db',varName,fYear[:4],result)
        write_sqlite_data(outdir+'/'+fYear+'.globalAveOcean.db',varName,fYear[:4],result)
    print(fYear,varDict,amocDict)   




atlmask = grid_info_from_gridspec("${gridspec}")
extract_ocean_scalar( str('${oname}'), '.','.', '')

EOF

#import sqlite3, cdms2, cdutil, MV2, numpy, cdtime
#import sys
## Set current year
#fYear = "${oname}"
## Test to see if sqlite databse exits, if not, then create it
#dbFile = "${localRoot}/db/.globalAveOcean.db"
## Read in ocean scalar annual file
#fdata = cdms2.open(fYear + '.ocean_scalar_month.nc')
#def extractScalarField(varName):
#  var=fdata(varName)
#  var=cdutil.YEAR(var).squeeze()
#  return var
##  return fdata(varName)[0]
#ignoreList = ['time_bounds', 'time_bnds', 'average_T2', 'average_T1', 'average_DT']
#varDict = fdata.variables
#varDict = list(set(varDict) - set(ignoreList))
#globalMeanDic={}
#for varName in varDict:
#  conn = sqlite3.connect("${localRoot}/db/.globalAveOcean.db")
#  c = conn.cursor()
#  globalMeanDic[varName] = extractScalarField(varName)
#  sql = 'create table if not exists ' + varName + ' (year integer primary key, value float)'
#  sqlres = c.execute(sql)
#  sql = 'insert or replace into ' + varName + ' values(' + fYear[:4] + ',' + str(globalMeanDic[varName]) + ')'
#  sqlres = c.execute(sql)
#  conn.commit()
#  c.close()
#  conn.close()
#EOF
cat > amoc.py <<EOF
import cdms2, netCDF4, tarfile, numpy, sqlite3
from scipy.io import netcdf
fYear = "${oname}"
#-- Attempt to read vmo and/or vh
doVmo = False
doVh = False
try:
  vmoFile = fYear+'.ocean_annual_rho2.nc'
  vmo = netCDF4.Dataset(vmoFile).variables['vmo'][0]
  vmo = vmo.filled(0)*1e-9
  zl  = netCDF4.Dataset(vmoFile).variables['rho2_l'][:]
  yq_vmo  = netCDF4.Dataset(vmoFile).variables['yq'][:]
  doVmo = True
except: doVmo = False
try:
  vhFile = fYear+'.ocean_z_month.nc'
  if 'vmo' in netCDF4.Dataset(vhFile).variables:
    vmo  = netCDF4.Dataset(vhFile).variables['vmo'][0]
    vmo  = vmo.filled(0)*1e-9
    zl  = netCDF4.Dataset(vhFile).variables['z_l'][:]
    yq_vmo  = netCDF4.Dataset(vhFile).variables['yq'][:]
    doVmo = True
  else:
    vh  = netCDF4.Dataset(vhFile).variables['vh'][0]
    vh  = vh.filled(0)*1e-9
    zt  = netCDF4.Dataset(vhFile).variables['zt'][:]
    yq_vh  = netCDF4.Dataset(vhFile).variables['yq'][:]
    doVh = True
except: doVh = False
#-- Get grid info from gridspec file
gsFile = "${gridspec}"
TF = tarfile.open(gsFile,'r')
member = [m for m in TF.getmembers() if 'ocean_hgrid' in m.name][0]
nc = netcdf.netcdf_file(TF.extractfile(member),'r')
x = nc.variables['x'][1::2,1::2]
y = nc.variables['y'][1::2,1::2]
member = [m for m in TF.getmembers() if 'topog' in m.name][0]
nc = netcdf.netcdf_file(TF.extractfile(member),'r')
depth = nc.variables['depth'][:]
def ice9it(i, j, depth, minD=0.):
  wetMask = 0*depth      
  (nj,ni) = wetMask.shape
  stack = set()
  stack.add( (j,i) )
  while stack:
    (j,i) = stack.pop()
    if wetMask[j,i] or depth[j,i] <= minD: continue
    wetMask[j,i] = 1
  
    if i>0: stack.add( (j,i-1) )
    else: stack.add( (j,ni-1) ) # Periodic beyond i=0
    if i<ni-1: stack.add( (j,i+1) )
    else: stack.add( (j,0) ) # Periodic beyond i=ni-1
    
    if j>0: stack.add((j-1,i))
        
    if j<nj-1: stack.add( (j+1,i) )
    else: stack.add( (j,ni-1-i) ) # Tri-polar fold beyond j=nj-1
  return wetMask
    
def ice9(x, y, depth, xy0):
  #ji = nearestJI(x, y, xy0)
  ji = nearestJI(x, y, xy0[0], xy0[1] )
  return ice9it(ji[1], ji[0], depth)
def nearestJI(x, y, x0, y0):
  return numpy.unravel_index( ((x-x0)**2 + (y-y0)**2).argmin() , x.shape)
def southOf(x, y, xy0, xy1):
  x0 = xy0[0]; y0 = xy0[1]; x1 = xy1[0]; y1 = xy1[1]
  dx = x1 - x0; dy = y1 - y0
  Y = (x-x0)*dy - (y-y0)*dx
  Y[Y>=0] = 1; Y[Y<=0] = 0
  return Y
#print 'Generating global wet mask ...',
wet = ice9(x, y, depth, (0,-35)) # All ocean points seeded from South Atlantic
#print 'done.'
code = 0*wet
#print 'Finding Cape of Good Hope ...',
tmp = 1 - wet; tmp[x<-30] = 0
tmp = ice9(x, y, tmp, (20,-30.))
yCGH = (tmp*y).min()
#print 'done.', yCGH
#print 'Finding Melbourne ...',
tmp = 1 - wet; tmp[x>-180] = 0
tmp = ice9(x, y, tmp, (-220,-25.))
yMel = (tmp*y).min()
#print 'done.', yMel
#print 'Processing Persian Gulf ...'
tmp = wet*( 1-southOf(x, y, (55.,23.), (56.5,27.)) )
tmp = ice9(x, y, tmp, (53.,25.))
code[tmp>0] = 11
wet = wet - tmp # Removed named points
#print 'Processing Red Sea ...'
tmp = wet*( 1-southOf(x, y, (40.,11.), (45.,13.)) )
tmp = ice9(x, y, tmp, (40.,18.))
code[tmp>0] = 10
wet = wet - tmp # Removed named points
#print 'Processing Black Sea ...'
tmp = wet*( 1-southOf(x, y, (26.,42.), (32.,40.)) )
tmp = ice9(x, y, tmp, (32.,43.))
code[tmp>0] = 7
wet = wet - tmp # Removed named points
#print 'Processing Mediterranean ...'
tmp = wet*( southOf(x, y, (-5.7,35.5), (-5.7,36.5)) )
tmp = ice9(x, y, tmp, (4.,38.))
code[tmp>0] = 6
wet = wet - tmp # Removed named points
#print 'Processing Baltic ...'
tmp = wet*( southOf(x, y, (8.6,56.), (8.6,60.)) )
tmp = ice9(x, y, tmp, (10.,58.))
code[tmp>0] = 9
wet = wet - tmp # Removed named points
#print 'Processing Hudson Bay ...'
tmp = wet*( 
           ( 1-(1-southOf(x, y, (-95.,66.), (-83.5,67.5)))
              *(1-southOf(x, y, (-83.5,67.5), (-84.,71.))) 
           )*( 1-southOf(x, y, (-70.,58.), (-70.,65.)) ) )
tmp = ice9(x, y, tmp, (-85.,60.))
code[tmp>0] = 8
wet = wet - tmp # Removed named points
#print 'Processing Arctic ...'
tmp = wet*( 
          (1-southOf(x, y, (-171.,66.), (-166.,65.5))) * (1-southOf(x, y, (-64.,66.4), (-50.,68.5))) # Lab Sea
     +    southOf(x, y, (-50.,0.), (-50.,90.)) * (1- southOf(x, y, (0.,65.5), (360.,65.5))  ) # Denmark Strait
     +    southOf(x, y, (-18.,0.), (-18.,65.)) * (1- southOf(x, y, (0.,64.9), (360.,64.9))  ) # Iceland-Sweden
     +    southOf(x, y, (20.,0.), (20.,90.)) # Barents Sea
     +    (1-southOf(x, y, (-280.,55.), (-200.,65.)))
          )
tmp = ice9(x, y, tmp, (0.,85.))
code[tmp>0] = 4
wet = wet - tmp # Removed named points
#print 'Processing Pacific ...'
tmp = wet*( (1-southOf(x, y, (0.,yMel), (360.,yMel)))
           -southOf(x, y, (-257,1), (-257,0))*southOf(x, y, (0,3), (1,3))
           -southOf(x, y, (-254.25,1), (-254.25,0))*southOf(x, y, (0,-5), (1,-5)) 
           -southOf(x, y, (-243.7,1), (-243.7,0))*southOf(x, y, (0,-8.4), (1,-8.4)) 
           -southOf(x, y, (-234.5,1), (-234.5,0))*southOf(x, y, (0,-8.9), (1,-8.9)) 
          )
tmp = ice9(x, y, tmp, (-150.,0.))
code[tmp>0] = 3
wet = wet - tmp # Removed named points
#print 'Processing Atlantic ...'
tmp = wet*(1-southOf(x, y, (0.,yCGH), (360.,yCGH)))
tmp = ice9(x, y, tmp, (-20.,0.))
code[tmp>0] = 2
wet = wet - tmp # Removed named points
#print 'Processing Indian ...'
tmp = wet*(1-southOf(x, y, (0.,yCGH), (360.,yCGH)))
tmp = ice9(x, y, tmp, (55.,0.))
code[tmp>0] = 5
wet = wet - tmp # Removed named points
#print 'Processing Southern Ocean ...'
tmp = ice9(x, y, wet, (0.,-55.))
code[tmp>0] = 1
wet = wet - tmp # Removed named points
code[wet>0] = -9
(j,i) = numpy.unravel_index( wet.argmax(), x.shape)
if j:
  print('There are leftover points unassigned to a basin code')
  while j:
    print (x[j,i],y[j,i],[j,i])
    wet[j,i]=0
    (j,i) = numpy.unravel_index( wet.argmax(), x.shape)
else: print('All points assigned a basin code')
#-- Define atlantic/arctic mask
atlmask = numpy.where(numpy.logical_or(code==2,code==4),1.,0.)

def MOCpsi(vh, vmsk=None, yq=None):
  shapeMask = list(vmsk.shape);
  shape = list(vh.shape); shape[-3] += 1
  if len(shape)==3:
    if vmsk is not None:
      if shape[1]==shapeMask[0]+1:
        # Symmetric grid (ignore southern most latitude)
        #print "vh_=vh[:,1:,:] shape1"
        vh_=vh[:,1:,:]
        yq_=yq[1:]
      elif shape[2]==shapeMask[1]+1:  
        vh_=vh[:,:,1:]
        yq_=yq[1:]
        #print "vh_=vh[:,:,1:] shape2"
      else:
        vh_ = vh
        yq_=yq
    shapevh = list(vh_.shape); 
    #print "shape vh_",shapevh
    #shapevh[-3] +=1
    psi = numpy.zeros(shapevh[:-1])
    for k in range(shapevh[-3]-1,0,-1):
      #print 'k=',k
      if vmsk is not None: 
        psi[k-1,:] = psi[k,:] - vh_[k-1].sum(axis=-1)
      else: 
        psi[k-1,:] = psi[k,:] - (vmsk*vh_[k-1]).sum(axis=-1)
  else:
    if vmsk is not None:
      if shape[2]==shapeMask[0]+1:  
        vh_=vh[:,:,1:,:]
        yq_=yq[1:]
      elif shape[3]==shapeMask[1]+1:  
        vh_=vh[:,:,:,1:]
        yq_=yq[1:]
      else:
        vh_=vh
    for n in range(shape[0]):
      for k in range(shape[-3]-1,0,-1):
        if vmsk==None: psi[n,k-1,:] = psi[n,k,:] - vh_[n,k-1].sum(axis=-1)
        else: psi[n,k-1,:] = psi[n,k,:] - (vmsk*vh_[n,k-1]).sum(axis=-1)
  return psi,yq_
varDict = {}
if doVmo == True:
  psi,yq_vmo_ = MOCpsi(vmo,vmsk=atlmask,yq=yq_vmo)
  print  ("shape of zl",list(zl.shape))
  print  ("shape of MOCpsi",list(psi.shape))
  print  ("shape of yq_vmo",list(yq_vmo.shape))
  print  ("shape of yq_vmo_",list(yq_vmo_.shape))
  maxsfn = numpy.max(psi[numpy.logical_and(zl>500,zl<2500)][:,numpy.greater_equal(yq_vmo_,20)])
  print('AMOC vmo = %s' % maxsfn)
  varDict['amoc_vmo'] = maxsfn
if doVh == True:
  psi,yq_vh_ = MOCpsi(vh,vmsk=atlmask,yq=yq_vh)
  maxsfn = numpy.max(psi[numpy.logical_and(zt>500,zt<2500)][:,numpy.greater_equal(yq_vh_,20)])
  print('AMOC vh = %s' % maxsfn)
  varDict['amoc_vh'] = maxsfn


# Comment out for now. 
#for varName in varDict:
#  conn = sqlite3.connect("${localRoot}/db/.globalAveOcean.db")
#  c = conn.cursor()
#  sql = 'create table if not exists ' + varName + ' (year integer primary key, value float)'
#  sqlres = c.execute(sql)
#  sql = 'insert or replace into ' + varName + ' values(' + fYear[:4] + ',' + str(varDict[varName]) + ')'
#  sqlres = c.execute(sql)
#  conn.commit()
#  c.close()
#  conn.close()
EOF

cat > merge.py << EOF
import sqlite3, sys
#
## usage python merge.py src dst
#
con = sqlite3.connect(sys.argv[2])
cur = con.cursor()
sql = "ATTACH '"+sys.argv[1]+"' as src"
cur.execute(sql)
cur.close()
cur = con.cursor()
sql = "SELECT * FROM main.sqlite_master WHERE type='table'"
cur.execute(sql)
main_tables = cur.fetchall()
cur.close()
cur = con.cursor()
sql = "SELECT * FROM src.sqlite_master WHERE type='table'"
cur.execute(sql)
src_tables = cur.fetchall()
cur.close()
for var in src_tables:
  varname = var[1]
  if varname not in [x[1] for x in src_tables]:
    cur = con.cursor()
    cur.execute(var[-1])
    cur.close()
  cur = con.cursor()
  sql = "INSERT OR REPLACE into "+varname+" SELECT * FROM src."+varname
  cur.execute(sql)
  cur.close()
con.commit()
con.close()
exit(0)
EOF
#-- Run the averager script
#gcp global*py gfdl:/home/$USER/tmp_refine
python global_atmos_ave.py
#python global_atmos_aer_ave.py
python global_ocean_ave.py
#python global_land_ave.py
#python amoc.py

ls -l *db

foreach reg (global nh sh tropics)
  foreach component (Atmos AtmosAer Land Ocean)
    if ( $component == "Ocean" && $reg != "global") then
      #only do OceanAveglobal
      break
    endif
    #if ($component == "Ice"   && $reg != "global") then
    #  #only do IceAveglobal
    #  break
    #endif

    if ( ! -f ${localRoot}/db/${reg}Ave${component}.db ) then
        echo "Copying " ${oname} ${reg}Ave${component}
        cp -fv ${oname}.${reg}Ave${component}.db ${localRoot}/db/${reg}Ave${component}.db
    else
        echo "Merging " ${oname} ${reg}Ave${component}
        python merge.py ${oname}.${reg}Ave${component}.db ${localRoot}/db/${reg}Ave${component}.db
    endif
  end 
end




echo "  ---------- end refineDiag_data_stager.csh ----------  "
echo ""
echo ""
echo ""
exit

