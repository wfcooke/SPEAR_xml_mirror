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


global_average_cubesphere(str('${oname}'), '.', '.', "Atmos", "atmos_month", '')
global_average_cubesphere(str('${oname}'), '.', '.', "AtmosAer", "atmos_month_aer", '')

global_average_land( str('${oname}'), '.', ".", "Land", "land_month", '')

#Move to global_oean_ave.py
