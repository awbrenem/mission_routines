




function BW_cntr_right,freqs,power,subscr,plotpow=plotpow,mVm=mVm,res=res
  ;Takes a plot of frequency vs. power (i.e. from an FFT) and calculates
  ;the bandwidth (BW) as double the distance from the frequency at
  ;maximum power to the highest frequency above the power cutoff.
  ;Note that I'm not counting any power at freqs lower than the peak freq
  ;because there is often a peak at low freqs that causes the bandwidth to be huge.

  e = 2.71828d0

  ;the cut off power is 1/e of the maximum power
  minpow = power[subscr]*(1.-1./e) ;power is on a linear scale
  minpow = power[subscr]*(1.-1./2.) ;power is on a linear scale



  ;save the data of frequency >= the frequency of maximum power
  power_tmp = power[subscr:n_elements(power)-1]  ;only data at and to the right of the peak
  ;Indices where the power is above the cut off power.
  tmp = where(power_tmp ge minpow) + subscr

  ;If anywhere is at least half the power of the maximum power, save it.
  if tmp[0] ne -1 then begin
    ;+
    ;bandwidth is twice the difference between the maximum frequency above the cut off
    ;and the minimum frequency below the cut off
    ;-
    BW = (max(freqs[tmp]) - min(freqs[tmp]))*2

    ;+
    ;Sometimes the bandwidth is zero b/c the power drops off to below the 1/e
    ;value within the frequency resolution. In this case,
    ;bandwidth = frequency resolution.
    ;-

    if BW eq 0 then BW = res
  endif else BW = !values.f_nan

  ;option to plot the power vs. frequency and which data points are above the cutoff
  if keyword_Set(plotpow) then begin
    plot,freqs,power-power[0],xrange=[0,5*freqs[subscr]]
    oplot,freqs[tmp],power[tmp]-power[0],psym=5,color=250
    wait,.5
  endif

  return,BW
END

function psp_narrowband_analyze, dates, time, wave, fce_data, type=type, fft_step=fft_step, step_overlap=step_overlap
  
  
  if ~keyword_set(fft_step) then fft_step = 0.1
  if ~keyword_set(step_overlap) then step_overlap = 0.5
  ;type = spec OR wavelet OR fft
  if ~keyword_set(type) then type = 'spec'
  
  notes = ["time_bins -> the time bins returned from slidespec or wavelet (secs)",$
    "freq_bins -> the freq bins returned from slidespec or wavelet (Hz)",$
    "max_power -> array of the maximum power for each time column",$
    "bandwidth -> array of bandwidth (Hz) as a function of time",$
    "timewidth -> array of the time of max power for each frequency row",$
    "freq_of_max_power -> array of the frequency of the max power for each time",$
    "df_f -> the bandwidth divided by the frequency for each time"]

  notes2 = ["max_power -> max power for the entire burst capture",$
    "bandwidth -> bandwidth for the time of maximum power in burst capture",$
    "freq_of_max_power -> the freq of the max power for the entire TDS capture",$
    "df_f -> the bandwidth divided by the frequency for the time of maximum power"]
  
  
  
  elems = n_elements(wave)
  dt = (max(time,/nan)-min(time,/nan))/n_elements(time)
  
  
  CASE type OF
    'spec':BEGIN

      ;sliding autocorrelation function that helps us id waves of high coherence
      ac = slide_autocor(time,wave,fft_step,step_overlap,time_index=time_index,amp_array=amp)
     
      power = slide_spec(time,wave,fft_step,step_overlap,iwindow=1,time_index=time_index,freq_bins=freqs)
      ;convert from kHz to Hz. Unnecessary unless you're commenting out the definition below
      
      freqs=freqs*1000.
      
      nfreqs = n_elements(freqs)
      ntimes = n_elements(time_index)
      namp = n_elements(amp)
      ;print, amp
      ;this is an array of the time elapsed at each step
      ;e.g. if time = [2.0,2.4,2.8,3.2,3.6,4.0], then telapsed = [0,0.4,0.8,1.2,1.6,2.0]
      telapsed = (max(time,/nan)-min(time,/nan))*indgen(n_elements(time_index))/(n_elements(time_index)-1.)

      ;empty struct to hold the data from a tds
      x = {time:dblarr(elems),$       ;array of times when samples are taken
        wave:dblarr(elems),$       ;values at each time
        time_bins:dblarr(ntimes),$     ;time since start at each index
        freq_bins:dblarr(nfreqs),$     ;frequency bins
        bandwidth:dblarr(ntimes-4),$     ;bandwidth of the max power freq at each step
        max_freq:dblarr(ntimes-4),$ ;maximum frequency in each bin
        power_maxima:dblarr(ntimes-4),$  ;max power in each bin
        df_f:dblarr(ntimes-4),$      ;normalized bandwidth at each time
        amp:dblarr(ntimes-4),$       ;working definition of amplitude (mV/m pk-pk)
        quality:dblarr(ntimes-4),$   ;working with Aarons definition of quality -------> (100*autocorrelate/df_f)*amplitude
        timelength:0.,$ ;length of signal above cutoff quality
        freq_maxima:0.,$ ;average frequency over length of the signal
        max_amp:0.,$ ;maximum amplitude of signal
        f_fce:0.,$ ;ratio of freq/fce
        bw_maxima:0.,$
        df_f_maxima:0.,$
        notes:notes $
      }
      
      
      
      x.time = time
      x.wave = wave
      x.time_bins = telapsed  ;s
      x.freq_bins = freqs   ;Hz
      ;help, power
      ;print, freqs[1:60]
      ;s = surface(power[*,1:60])
      ;print, ' '
      ;print, fce_data*0.15
      
      a = where((freqs lt (1.6)*fce_data) and (freqs gt 20))
      
      ;help, max(power[*,a])
      ;help, a, freqs, fce_data
      ;print,a, freqs[50:500], fce_data[0]/2
      for tt=2,ntimes-3 do begin
        x.power_maxima[tt-2] = max(power[tt,a], sub)
        x.max_freq[tt-2] = freqs[a[sub]]
        x.bandwidth[tt-2] = BW_cntr_right(freqs[a],power[tt,a],sub,res=(freqs[1]-freqs[0]))
        x.df_f[tt-2] = x.bandwidth[tt-2]/x.max_freq[tt-2]
        x.amp[tt-2] = amp[tt]
        x.quality[tt-2] = (100*ac[tt]/x.bandwidth[tt-2])*power[tt,a[sub]]
      endfor
      ;3p = plot(x.quality)
      
      b = where((x.quality gt 0.01))
      startpoint = b[0]
      endpoint = b[n_elements(b)-1]

      duration = abs(x.time_bins[startpoint+2]-x.time_bins[endpoint+2])

      res = telapsed[1]

      if b[0] eq -1 then begin
        x.timelength = -1
      endif else begin
        if duration eq 0 then x.timelength = res $
        else x.timelength = duration + res
      endelse 
      
      max_power = max(power[2:ntimes-3,a], subscr2)
      
      indices = ARRAY_INDICES(power[2:ntimes-3,a], subscr2)
      subtime = indices[0]
      
      if n_elements(indices) lt 2 then begin
        x.freq_maxima = 0
        x.quality[*] = -1
        print, 'Index error: skipping this wave'
      endif else begin
        subfreq = a[indices[1]]
        x.freq_maxima = freqs[subfreq]
      endelse
     
      x.bw_maxima = x.bandwidth[subtime]
      ;if x.bw_maxima gt 11 then x.quality[*]=-1      
      x.df_f_maxima = x.df_f[subtime]     
      x.f_fce = x.freq_maxima/fce_data     
      x.max_amp = max(x.amp); max(x.amp[b])     
      
      for i=0,n_elements(x.quality)-1 do begin
        ;-------------------;
        ;   Test for Dust   ;
        ;-------------------;

        dust = 'no'

        goo1 = where(wave gt 0.)
        goo2 = where(wave lt 0.)

        if goo1[0] ne -1 then max1 = max(wave[goo1]) else max1 = 0.
        if goo2[0] ne -1 then max2 = min(wave[goo2]) else max2 = 0.

        if max1 ne 0. and max2 ne 0. then begin
          if abs(max1) ge abs(max2) then wr = abs(max1/max2) else wr = abs(max2/max1)
        endif else wr = 10. ;arbitrary number to classify as dust


        t1 = where(ac[i+2] gt 0.)
        t2 = where(ac[i+2] lt 0.)

        if t1[0] ne -1 and t2[0] ne -1 then begin
          max1 = max(ac(t1))
          max2 = min(ac(t2))
          if abs(max1) ge abs(max2) then ar = abs(max1/max2) else ar = abs(max2/max1)
        endif else begin
          ar = 1000.  ;arbitrary ratio
        endelse

        ;four conditions that classify wave as dust
        ;if any quality is less than 1, it's dust
        ;THIS NEEDS TO BE REDEFINED 10/20
        check = boolean(1)
        ;for zz=0,3 do begin

        check = check && (x.quality[0] le 1.)

        ;endfor

        ;if check then dust = 'yes'
        s = where(x.quality[i] gt 20000)
        if s[0,0] ne -1 then dust = 'yes'
        ;print, dust, ' 1'
;        if wr ge 7 and ar ge 3 then dust = 'yes'
;        print, dust, ' 2'
;        if ar ge 7 and wr ge 3 then dust = 'yes'
;        print,dust, ' 3'
;        if ar eq 1000 then dust = 'yes'
;        print, dust, ' 4'
        ;if x.max_amp gt 75 and x.timelength lt 0.1 then dust = 'yes'
        ;if x.max_amp gt 180 and x.timelength lt 0.2 then dust = 'yes'
        if x.max_amp gt 1500 then dust = 'yes' 
        
        ;print, dust, ' 2'
        ;On the last pass through, compare amplitudes in separate channels as well as
        ;absolute amplitude.
        if (max(abs(wave)) gt 400) then dust='yes'
        
        ;print, dust, ' 3'
        
        if dust eq 'yes' then x.quality[*] = -1

      endfor
      ;print, x.quality
      
  ENDCASE
end

return, x

end