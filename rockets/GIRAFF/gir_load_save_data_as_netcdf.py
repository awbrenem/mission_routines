"""
Save GIRAFF data as netcdf. This can be loaded in IDL

We can then load the data in IDL with:
path = '/Users/abrenema/Desktop/Research/Rocket_missions/GIRAFF/data/efield_VLF/'
fn = 'Giraff_381_efield_vlf_calibrated.netcdf'

ncid = NCDF_OPEN(path+fn)
E12id = NCDF_VARID(ncid, 'E12')
E34id = NCDF_VARID(ncid, 'E34')
timeid = NCDF_VARID(ncid, 'time')

;get a variable out of the data
NCDF_VARGET, ncid, E12id, E12
NCDF_VARGET, ncid, E34id, E34
NCDF_VARGET, ncid, timeid, time

;Close the netcdf file
NCDF_CLOSE, ncid


time = time + time_double('1980-01-01/00:00')
store_data,'E12',time,E12
store_data,'E34',time,E34

"""


import sys 
sys.path.append('/Users/abrenema/Desktop/code/Aaron/github/mission_routines/rockets/GIRAFF/')
from gir_load_fields import GIRAFF_Fields_Loader as GFL
import netCDF4
import numpy as np


pld = '381'
#v12 = GFL(pld,'VLF12D')
#v34 = GFL(pld,'VLF34D')
v12 = GFL(pld,'V12D')
v34 = GFL(pld,'V34D')

fs = v12.chnspecs['fs']
wf12, tvals = v12.load_data()
wf34, tgoo = v34.load_data()
#wf12, tvals = v12.load_data_gainphase_corrected()
#wf34, tgoo = v34.load_data_gainphase_corrected()



#--------------------------------------------
#Save data to a Netcdf file
#see tutorial at https://unidata.github.io/python-training/workshop/Bonus/netcdf-writing/
#--------------------------------------------

save_dir = '/Users/abrenema/Desktop/'
filename = 'datatst.netcdf'


nc = netCDF4.Dataset(save_dir + filename,'w',format='NETCDF4')


time_dim = nc.createDimension('time',len(tvals))
E12_dim = nc.createDimension('E12',len(wf12))
E34_dim = nc.createDimension('E34',len(wf12))


#independent variables
time = nc.createVariable('time',np.double, ('time'))
time.units = 'sec'
time.long_name = 'sec since launch'



#dependent variables
E12 = nc.createVariable('E12',np.float32, ('time'))
E12.standard_name = 'E12 mVm'
E12.long_name = 'GIRAFF calibrated E12 (mV/m)'
E34 = nc.createVariable('E34',np.float32, ('time'))
E34.standard_name = 'E34 mVm'
E34.long_name = 'GIRAFF calibrated E34 (mV/m)'


#assign data
E12[:] = wf12
E34[:] = wf34
time[:] = tvals

nc.close()

print('h')