

pro download_to_extdrive
   
   if getenv('PSP_STAGING_ID') EQ '' then setenv, 'PSP_STAGING_ID=short186'
   if getenv('PSP_STAGING_PW') EQ '' then setenv, 'PSP_STAGING_PW=poLuSHISheoc'
   
   get_timespan, time
   ;                       
   setenv, 'ROOT_DATA_DIR=/Volumes/500GB/Users/benshort/data/'
   
   spp_fld_load, trange=time, type='dfb_dbm_dvac',/no_load           
   
end
