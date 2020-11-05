;+
;FUNCTION:   slide_autocor.pro
;  Calculates a sliding autocorrelation for time series data
;
;
;  version 0.1, 20190912
;
;ARGUMENTS:
;   TIME    -> time array (sec)
;   DATA    -> time series data array
;   STEP    -> length of each slice, in percent of total length
;   STEP_OVERLAP  -> slice overlap, given in percent of slice length
;
;KEYWORDS:
;       TIME_INDEX=time_index <- integer array containing the index for each slice time.
;                  i.e., time[time_index] is the array of x-axis times for the spec
;   
;   
;               if not set then: power_temp=2*padded_step_length*(abs(a_correlate(temp_time_series, templag)))^2
;                 
;   iwindow         -> window the data (1:Hamming, 2:Hanning, 3:Gaussian)
;   zero_pad        -> set to have program zero pad the data. Can run MUCH faster if you set this.
;
;RETURNS:
;   AC_ARRAY <- 1D sliding autocorrelate array (# time steps)
;
;CALLING SEQUENCE:
;       IDL>
;
;NOTES:
;       (none)
;-
;CREATED BY:    Ben Short, 9-12-2019
;
;HISTORY:
;   
;
;
;INCLUDED MODULES:
;       
;
;LIBS USED:
;       (none)
;
;DEPENDENCIES
;
;-

function slide_autocor,time,time_series,step,step_overlap,time_index=time_index,amp_array=amp_array

  npoints=0L    ; number of points in the time series
  step_points=0L  ; number of points in each "slice"
  step_overlap_points=0L
  step_start=0L ; starting index of sliding "slice"
  step_end=0L   ; ending index of sliding "slice"
  step_count=0L ; counter for stepping through the time series

  ; figure out number of points per step
  if step gt 1. then step=1.
  if step_overlap gt 1. then step_overlap=1.
  if step le 0. then step=.25
  npoints=size(time,/n_elements)
  if npoints ne size(time_series,/n_elements) then message,'size mismatch in time series.'
  step_points=long(npoints*step)
  step_overlap_points=long(step_points*step_overlap)
  ; make sure we have a working overlap (i.e., in the range 0,step_points-1)
  if step_overlap_points lt 0 then step_overlap_points=0
  if step_overlap_points eq step_points then step_overlap_points=step_points-1

  ; figure out time step (assuming input times in seconds)
  total_time=time[npoints-1]-time[0]
  time_step=total_time/npoints
  step_length=step_points*time_step

  step_start=0
  step_end=step_points-1
  
  padded_step_points=0L
  
  padded_step_points=step_points
  padded_step_length=step_length
  
  
  freq_bins=(lindgen(padded_step_points/2-1)+1)*npoints/(padded_step_points*total_time*1000.)
  
  window_array=make_array(step_points,value=1.0,/float)
  
  sumsq=total(window_array^2)

  rms=sqrt(sumsq/step_points)
  
  window_array=window_array/rms
  
  
  ac_array=fltarr((npoints-step_overlap_points)/(step_points-step_overlap_points))
  amp_array=fltarr((npoints-step_overlap_points)/(step_points-step_overlap_points))
  
  while step_end le npoints-1 do begin
    
    temp_time_series=window_array*time_series[step_start:step_end]
    
    temp_lag = 2*indgen(n_elements(temp_time_series),/LONG) - (n_elements(temp_time_series) - 1)
    
    ac_temp=a_correlate(temp_time_series/max(temp_time_series, /nan), temp_lag,/double)
    
    half = n_elements(ac_temp)/2.
    actmp = reverse(ac_temp[0:half-1])

    jj = 0l
    maxim = 0l

    gt0 = where(actmp gt 0.)
    
    if gt0[0] ne -1 then begin
      ;find first zero crossing of autocorrelation
      while (actmp[jj] ge 0.) and (jj lt n_elements(actmp)-2) do jj++
      ;get back to positive autocorrelation values
      while (actmp[jj] lt maxim) and (jj lt n_elements(actmp)-2) do jj++
      ;find max positive peak from rest of autocorrelation array
    endif
    
    maxim = max(actmp[jj:n_elements(actmp)-1],sub)

    ratio_ac = abs(maxim/actmp[0])  ;ratio of max peak (value=1 for autocorrelation) to the second peak
    
    amp_temp=abs(max(temp_time_series)-min(temp_time_series))

    ac_array[step_count]=ratio_ac
    amp_array[step_count] = amp_temp
    step_count+=1
    step_start+=(step_points-step_overlap_points)
    step_end+=(step_points-step_overlap_points)
  endwhile
  
  ;amplitude array
  ;amp_array=fltarr((npoints-step_overlap_points)/(step_points-step_overlap_points))
  ;while step_end le npoints-1 do begin

   ; temp_time_series=window_array*time_series[step_start:step_end]
     
   ; amp_temp=max(temp_time_series)-min(temp_time_series)
    
    ;print, amp_temp
    
    ;amp_array[stepcount] = amp_temp
   ; step_count+=1
    ;step_start+=(step_points-step_overlap_points)
    ;step_end+=(step_points-step_overlap_points)

;  endwhile
  
  
  
  ; set up the time index with padding
  time_index=lindgen(step_count)*(step_points-step_overlap_points)+step_points/2
  start_pad_index=(lindgen(round(float(time_index[0])/(step_points-step_overlap_points))))*(step_points-step_overlap_points)
  end_pad_index=(lindgen((npoints-time_index[step_count-1])/(step_points-step_overlap_points)))*(step_points-step_overlap_points)+time_index[step_count-1]
  time_index=[start_pad_index,time_index[0],time_index,end_pad_index,npoints-1]
  
  ;make sure to skip the two time indices on each end.
  ; pad the array
  start_pad_data=make_array(size(start_pad_index,/n_elements)+1,value=0)
  end_pad_data=make_array(size(end_pad_index,/n_elements)+1,value=0)
  ac_array=[start_pad_data,ac_array,end_pad_data]
  amp_array=[start_pad_data,amp_array,end_pad_data]
  
  
  return, ac_array 
  
end
  