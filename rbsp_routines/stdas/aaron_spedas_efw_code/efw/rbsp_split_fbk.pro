;+
; NAME:	rbsp_split_fbk
;
;
; PURPOSE:	Split the filterbank data into separate tplot variables for
;			each channel. Can also combine peak and average on a single panel
;
; KEYWORDS:
;	probe = 'a' or 'b'  NOTE: single spacecraft only, does not accept ['a b']
;	combine -> set to combine peak and average onto a single plot
;	meansz -> number of data points over which to calculate the number of values
;			  the local mean is away from the local standard deviation. This is used
;			  to set y-scaling. If not done then the yscale of the FBK plots
;			  is often dominated by single large amplitude spikes, masking the majority
;			  of the data.  Default is 100.
;;	meansz -> number of data points over which to calculate the mean. This is used
;;			  to set y-scaling. If not done then the yscale of the FBK plots
;;			  is often dominated by single large amplitude spikes, masking the majority
;;			  of the data.  Default is 20.
;;	ysc -> Scale factor for y-scaling. Default is 1.  This should be set to
;;			values greater than 1 for meansz larger than the default.  For
;;			example, ysc=3. works well for meansz=1000.
;
; CREATED: Aaron Breneman 11/07/2012
;
; MODIFIED: changed y-scaling based on a running average rather than the max value
;			for each FBK channel. This avoids having a ridiculous yscaling based on a few
;			very large amplitude spiky events.
;
; VERSION:
;$LastChangedBy: aaronbreneman $
;$LastChangedDate: 2014-04-21 15:19:06 -0500 (Mon, 21 Apr 2014) $
;$LastChangedRevision: 14901 $
;$URL: svn+ssh://thmsvn@ambrosia.ssl.berkeley.edu/repos/spdsoft/trunk/general/missions/rbsp/efw/rbsp_split_fbk.pro $
;
;-



pro rbsp_split_fbk,probe,combine=combine,meansz=sz,verbose=verbose,ysc=ysc


	vb = keyword_set(verbose) ? verbose : 0
	vb = vb > !rbsp_efw.verbose
	start_time=systime(1)

	;Frequency bins for FBK13 and FBK7
	fbk13_bins=['0.8-1.5', '1.5-3', '3-6', '6-12', '12-25', '25-50', $
				'50-100', '100-200', '200-400', '400-800', $
				'800-1600', '1600-3200', '3200-6500']
	fbk7_bins=fbk13_bins[lindgen(7)*2]




	;Determine ylim based on the mean value. Doing this avoids large y-scalings
	;based on very spiky FBK data
	if ~keyword_set(sz) then sz = 100.


	;the ratio of the local value of the peak FBK value divided by the standard deviation
	if ~keyword_set(maxrat) then maxrat = 30


	;grab dlimits structure for each of the channels to determine the source
	get_data,'rbsp'+probe+'_efw_fbk_7_fb1_pk',dlimits=fb7_fb1
	get_data,'rbsp'+probe+'_efw_fbk_7_fb2_pk',dlimits=fb7_fb2
	get_data,'rbsp'+probe+'_efw_fbk_13_fb1_pk',dlimits=fb13_fb1
	get_data,'rbsp'+probe+'_efw_fbk_13_fb2_pk',dlimits=fb13_fb2


	;Split up FBK13-FB1 peak array, if it exists
	get_data,'rbsp'+probe +'_efw_fbk_13_fb1_pk',data=goo,dlimits=dlim
	tst = size(goo,/dimensions)
	if tst ne 0 then begin
		for j=0,12 do store_data,'rbsp'+probe +'_fbk1_13pk_'+strtrim(j,2),data={x:goo.x,y:goo.y[*,j]}
		units = dlim.data_att.units
		for j=0,12 do options,'rbsp'+probe +'_fbk1_13pk_'+strtrim(j,2),'ytitle','pk!C'+fbk13_bins[j]+'!CHz!C['+units+']'
		options,'rbsp'+probe +'_fbk1_13??_*','labels','RBSP'+probe+'!C  '+fb13_fb1.data_att.channel
	endif

	get_data,'rbsp'+probe +'_efw_fbk_13_fb2_pk',data=goo,dlimits=dlim
	tst = size(goo,/dimensions)
	if tst ne 0 then begin
		for j=0,12 do store_data,'rbsp'+probe +'_fbk2_13pk_'+strtrim(j,2),data={x:goo.x,y:goo.y[*,j]}
		units = dlim.data_att.units
		for j=0,12 do options,'rbsp'+probe +'_fbk2_13pk_'+strtrim(j,2),'ytitle','pk!C'+fbk13_bins[j]+'!CHz!C['+units+']'
		options,'rbsp'+probe +'_fbk2_13??_*','labels','RBSP'+probe+'!C  '+fb13_fb2.data_att.channel
	endif

	get_data,'rbsp'+probe +'_efw_fbk_7_fb1_pk',data=goo,dlimits=dlim
	tst = size(goo,/dimensions)
	if tst ne 0 then begin
		units = dlim.data_att.units
		for j=0,6 do store_data,'rbsp'+probe +'_fbk1_7pk_'+strtrim(j,2),data={x:goo.x,y:goo.y[*,j]}
		for j=0,6 do options,'rbsp'+probe +'_fbk1_7pk_'+strtrim(j,2),'ytitle','pk!C'+fbk7_bins[j]+'!CHz!C['+units+']'
		options,'rbsp'+probe +'_fbk1_7??_*','labels','RBSP'+probe+'!C  '+fb7_fb1.data_att.channel
	endif

	get_data,'rbsp'+probe +'_efw_fbk_7_fb2_pk',data=goo,dlimits=dlim
	tst = size(goo,/dimensions)
	if tst ne 0 then begin
		units = dlim.data_att.units
		for j=0,6 do store_data,'rbsp'+probe +'_fbk2_7pk_'+strtrim(j,2),data={x:goo.x,y:goo.y[*,j]}
		for j=0,6 do options,'rbsp'+probe +'_fbk2_7pk_'+strtrim(j,2),'ytitle','pk!C'+fbk7_bins[j]+'!CHz!C['+units+']'
		options,'rbsp'+probe +'_fbk2_7??_*','labels','RBSP'+probe+'!C  '+fb7_fb2.data_att.channel
	endif

	get_data,'rbsp'+probe +'_efw_fbk_13_fb1_av',data=goo,dlimits=dlim
	tst = size(goo,/dimensions)
	if tst ne 0 then begin
		units = dlim.data_att.units
		for j=0,12 do store_data,'rbsp'+probe +'_fbk1_13av_'+strtrim(j,2),data={x:goo.x,y:goo.y[*,j]}
		for j=0,12 do options,'rbsp'+probe +'_fbk1_13av_'+strtrim(j,2),'ytitle','av!C'+fbk13_bins[j]+'!CHz!C['+units+']'
		options,'rbsp'+probe +'_fbk1_13??_*','labels','RBSP'+probe+'!C  '+fb13_fb1.data_att.channel
	endif

	get_data,'rbsp'+probe +'_efw_fbk_13_fb2_av',data=goo,dlimits=dlim
	tst = size(goo,/dimensions)
	if tst ne 0 then begin
		units = dlim.data_att.units
		for j=0,12 do store_data,'rbsp'+probe +'_fbk2_13av_'+strtrim(j,2),data={x:goo.x,y:goo.y[*,j]}
		for j=0,12 do options,'rbsp'+probe +'_fbk2_13av_'+strtrim(j,2),'ytitle','av!C'+fbk13_bins[j]+'!CHz!C['+units+']'
		options,'rbsp'+probe +'_fbk2_13??_*','labels','RBSP'+probe+'!C  '+fb13_fb2.data_att.channel
	endif

	get_data,'rbsp'+probe +'_efw_fbk_7_fb1_av',data=goo,dlimits=dlim
	tst = size(goo,/dimensions)
	if tst ne 0 then begin
		units = dlim.data_att.units
		for j=0,6 do store_data,'rbsp'+probe +'_fbk1_7av_'+strtrim(j,2),data={x:goo.x,y:goo.y[*,j]}
		for j=0,6 do options,'rbsp'+probe +'_fbk1_7av_'+strtrim(j,2),'ytitle','av!C'+fbk7_bins[j]+'!CHz!C['+units+']'
		options,'rbsp'+probe +'_fbk1_7??_*','labels','RBSP'+probe+'!C  '+fb7_fb1.data_att.channel
	endif

	get_data,'rbsp'+probe +'_efw_fbk_7_fb2_av',data=goo,dlimits=dlim
	tst = size(goo,/dimensions)
	if tst ne 0 then begin
		units = dlim.data_att.units
		for j=0,6 do store_data,'rbsp'+probe +'_fbk2_7av_'+strtrim(j,2),data={x:goo.x,y:goo.y[*,j]}
		for j=0,6 do options,'rbsp'+probe +'_fbk2_7av_'+strtrim(j,2),'ytitle','av!C'+fbk7_bins[j]+'!CHz!C['+units+']'
		options,'rbsp'+probe +'_fbk2_7??_*','labels','RBSP'+probe+'!C  '+fb7_fb2.data_att.channel
	endif



	if keyword_set(combine) then begin

		;Combine FBK pk and av quantities onto a single plot
		tst = tnames('rbsp'+probe +'_fbk1_13pk_0')
		if tst ne '' then begin
			for j=0,12 do store_data,'rbsp'+probe +'_fbk1_13comb_'+strtrim(j,2),data=['rbsp'+probe +'_fbk1_13pk_'+strtrim(j,2),'rbsp'+probe +'_fbk1_13av_'+strtrim(j,2)]
			options,'rbsp'+probe +'_fbk1_13comb_*','labels','RBSP'+probe+'!C  '+fb13_fb1.data_att.channel
		endif

		tst = tnames('rbsp'+probe +'_fbk2_13pk_0')
		if tst ne '' then begin
			for j=0,12 do store_data,'rbsp'+probe +'_fbk2_13comb_'+strtrim(j,2),data=['rbsp'+probe +'_fbk2_13pk_'+strtrim(j,2),'rbsp'+probe +'_fbk2_13av_'+strtrim(j,2)]
			options,'rbsp'+probe +'_fbk2_13comb_*','labels','RBSP'+probe+'!C  '+fb13_fb2.data_att.channel
		endif

		tst = tnames('rbsp'+probe +'_fbk1_7pk_0')
		if tst ne '' then begin
			for j=0,6 do store_data,'rbsp'+probe +'_fbk1_7comb_'+strtrim(j,2),data=['rbsp'+probe +'_fbk1_7pk_'+strtrim(j,2),'rbsp'+probe +'_fbk1_7av_'+strtrim(j,2)]
			options,'rbsp'+probe +'_fbk1_7comb_*','labels','RBSP'+probe+'!C  '+fb7_fb1.data_att.channel
		endif

		tst = tnames('rbsp'+probe +'_fbk2_7pk_0')
		if tst ne '' then begin
			for j=0,6 do store_data,'rbsp'+probe +'_fbk2_7comb_'+strtrim(j,2),data=['rbsp'+probe +'_fbk2_7pk_'+strtrim(j,2),'rbsp'+probe +'_fbk2_7av_'+strtrim(j,2)]
			options,'rbsp'+probe +'_fbk2_7comb_*','labels','RBSP'+probe+'!C  '+fb7_fb2.data_att.channel
		endif

		;Change colors to black and red
		tst = tnames('rbsp'+probe +'_fbk1_13comb_0')
		if tst ne '' then options,['rbsp'+probe+'_fbk1_13comb_?','rbsp'+probe+'_fbk1_13comb_??'],colors=[0,6]

		tst = tnames('rbsp'+probe +'_fbk2_13comb_0')
		if tst ne '' then options,['rbsp'+probe+'_fbk2_13comb_?','rbsp'+probe+'_fbk2_13comb_??'],colors=[0,6]

		tst = tnames('rbsp'+probe +'_fbk1_7comb_0')
		if tst ne '' then options,'rbsp'+probe+'_fbk1_7comb_?',colors=[0,6]

		tst = tnames('rbsp'+probe +'_fbk2_7comb_0')
		if tst ne '' then options,'rbsp'+probe+'_fbk2_7comb_?',colors=[0,6]


		tst = tnames('rbsp'+probe +'_fbk1_13comb_0')
		get_data,'rbsp'+probe+'_efw_fbk_13_fb1_pk',dlimits=dlim,data=goo
		if tst ne '' then begin
			units = dlim.data_att.units
			for j=0,12 do options,'rbsp'+probe +'_fbk1_13comb_'+strtrim(j,2),'ytitle',fbk13_bins[j]+'!CHz!C['+units+']'
		endif

		tst = tnames('rbsp'+probe +'_fbk2_13comb_0')
		get_data,'rbsp'+probe+'_efw_fbk_13_fb2_pk',dlimits=dlim,data=goo
		if tst ne '' then begin
			units = dlim.data_att.units
			for j=0,12 do options,'rbsp'+probe +'_fbk2_13comb_'+strtrim(j,2),'ytitle',fbk13_bins[j]+'!CHz!C['+units+']'
		endif

		tst = tnames('rbsp'+probe +'_fbk1_7comb_0')
		get_data,'rbsp'+probe+'_efw_fbk_7_fb1_pk',dlimits=dlim,data=goo
		if tst ne '' then begin
			units = dlim.data_att.units
			for j=0,6 do options,'rbsp'+probe +'_fbk1_7comb_'+strtrim(j,2),'ytitle',fbk7_bins[j]+'!CHz!C['+units+']'
		endif

		tst = tnames('rbsp'+probe +'_fbk2_7comb_0')
		get_data,'rbsp'+probe+'_efw_fbk_7_fb2_pk',dlimits=dlim,data=goo
		if tst ne '' then begin
			units = dlim.data_att.units
			for j=0,6 do options,'rbsp'+probe +'_fbk2_7comb_'+strtrim(j,2),'ytitle',fbk7_bins[j]+'!CHz!C['+units+']'
		endif
	endif


;The number of y-ticks is usually too large when 7 or 13 plots are displayed at once.
;Unfortunately reducing the tick numbers means that negative numbers show up. So, I have
;to manually set the yrange
options,'rbsp'+probe+'_fbk?_13pk_*','yticks',2
options,'rbsp'+probe+'_fbk?_7pk_*','yticks',2

;If yticks are manually set, yrange isn't set or is set to a range of zero (max=min),
;ystyle is set to 1, and the data is a horizontal line (max(y)=min(y)), then the IDL
;PLOT routine crashes. This next bit checks to see if any of the pk vars are horizontal
;lines, and if so, sets yticks to 0 (automatic) for that variable.
;--------------
pk13vars=tnames('rbsp'+probe+'_fbk1_13pk_*',cnt)
if cnt gt 0 then begin
	for i=0, n_elements(pk13vars)-1 do begin
		var=pk13vars[i]
		get_data,var,data=vardat
		if min(vardat.y) eq max(vardat.y) then $
		options,var,'yticks',0
	endfor
endif
pk13vars=tnames('rbsp'+probe+'_fbk2_13pk_*',cnt)
if cnt gt 0 then begin
	for i=0, n_elements(pk13vars)-1 do begin
		var=pk13vars[i]
		get_data,var,data=vardat
		if min(vardat.y) eq max(vardat.y) then $
		options,var,'yticks',0
	endfor
endif
pk7vars=tnames('rbsp'+probe+'_fbk1_7pk_*',cnt)
if cnt gt 0 then begin
	for i=0, n_elements(pk7vars)-1 do begin
		var=pk7vars[i]
		get_data,var,data=vardat
		if min(vardat.y) eq max(vardat.y) then $
		options,var,'yticks',0
	endfor
endif
pk7vars=tnames('rbsp'+probe+'_fbk2_7pk_*',cnt)
if cnt gt 0 then begin
	for i=0, n_elements(pk7vars)-1 do begin
		var=pk7vars[i]
		get_data,var,data=vardat
		if min(vardat.y) eq max(vardat.y) then $
		options,var,'yticks',0
	endfor
endif


;-----------------------------------------------
;Set ylimit and tick format based on keeping peak values within a certain number
;of standard deviations. Avoids having the scaling set by spiky noise
;-----------------------------------------------

;---------
;FBK1 - 13



tst = tnames('rbsp'+probe +'_fbk1_13pk_0')
if tst ne '' then begin
	get_data,'rbsp'+probe +'_fbk1_13pk_0',data=goo0

	nchunks = floor(n_elements(goo0.y)/sz)
	mom = fltarr(nchunks,13)
	maxv = fltarr(nchunks,13) ;peak values
	varv = fltarr(nchunks,13) ;variances

	ylim_start_time=systime(1)
	dmessage="using nchunks in fbk1_13pk: "+string(nchunks)
	dprint,verbose=verbose,dlevel=2,dmessage

	for j=0,12 do begin
		get_data,'rbsp'+probe +'_fbk1_13pk_'+strtrim(j,2),data=goob
		for i=0L,nchunks-1*sz do begin
			maxv[i,j] = max(goob.y[i*sz:i*sz+sz])
			tst = moment(goob.y[i*sz:i*sz+sz])
			varv[i,j] = sqrt(tst[1])
			mom[i,j]=total(goob.y[i*sz:i*sz+sz])/sz
		endfor
	endfor

	stdevv = sqrt(varv)
	ratio = maxv/stdevv


	;Reject values (for the y-scaling) where the max value is >> than the standard deviation
	good = bytarr(nchunks,13)
	good[*] = 1b

	for j=0,12 do begin
		goo = where(ratio[*,j] ge maxrat)
		if goo[0] ne -1 then good[goo,j] = 0b
	endfor

	good = float(good)

	for j=0,12 do ylim,'rbsp'+probe +'_fbk1_13pk_'+strtrim(j,2),0.,max(maxv[*,j]*good[*,j]),0
	if keyword_set(combine) then for j=0,12 do ylim,'rbsp'+probe +'_fbk1_13comb_'+strtrim(j,2),0,max(maxv[*,j]*good[*,j]),0

	dmessage='ylim fbk1_13pk runtime (s): '+string(systime(1)-ylim_start_time)
	dprint,verbose=verbose,dlevel=2,dmessage


	;Keep only 1 number after decimal place. More just clutters up plot
	for j=0,12 do options,'rbsp'+probe +'_fbk1_13pk_'+strtrim(j,2),'ytickformat','(f6.1)'
endif

;---------
;FBK2 - 13
;---------

tst = tnames('rbsp'+probe +'_fbk2_13pk_0')
if tst ne '' then begin
	get_data,'rbsp'+probe +'_fbk2_13pk_0',data=goo0

	nchunks = floor(n_elements(goo0.y)/sz)
	mom = fltarr(nchunks,13)
	maxv = fltarr(nchunks,13) ;peak values
	varv = fltarr(nchunks,13) ;variances

	ylim_start_time=systime(1)
	dmessage="using nchunks in fbk2_13pk: "+string(nchunks)
	dprint,verbose=verbose,dlevel=2,dmessage

	for j=0,12 do begin
		get_data,'rbsp'+probe +'_fbk2_13pk_'+strtrim(j,2),data=goob
		for i=0L,nchunks-1*sz do begin
			maxv[i,j] = max(goob.y[i*sz:i*sz+sz])
			tst = moment(goob.y[i*sz:i*sz+sz])
			varv[i,j] = sqrt(tst[1])
			mom[i,j]=total(goob.y[i*sz:i*sz+sz])/sz
		endfor
	endfor

	stdevv = sqrt(varv)
	ratio = maxv/stdevv


	;Reject values (for the y-scaling) where the max value is >> than the standard deviation
	good = bytarr(nchunks,13)
	good[*] = 1b

	for j=0,12 do begin
		goo = where(ratio[*,j] ge maxrat)
		if goo[0] ne -1 then good[goo,j] = 0b
	endfor

	good = float(good)

	for j=0,12 do ylim,'rbsp'+probe +'_fbk2_13pk_'+strtrim(j,2),0.,max(maxv[*,j]*good[*,j]),0
	if keyword_set(combine) then for j=0,12 do ylim,'rbsp'+probe +'_fbk2_13comb_'+strtrim(j,2),0,max(maxv[*,j]*good[*,j]),0

	dmessage='ylim fbk2_13pk runtime (s): '+string(systime(1)-ylim_start_time)
	dprint,verbose=verbose,dlevel=2,dmessage


	;Keep only 1 number after decimal place. More just clutters up plot
	for j=0,12 do options,'rbsp'+probe +'_fbk2_13pk_'+strtrim(j,2),'ytickformat','(f6.1)'
endif

;---------
;FBK1 - 7
;---------

tst = tnames('rbsp'+probe +'_fbk1_7pk_0')
if tst ne '' then begin
	get_data,'rbsp'+probe +'_fbk1_7pk_0',data=goo0

	nchunks = floor(n_elements(goo0.y)/sz)
	mom = fltarr(nchunks,7)
	maxv = fltarr(nchunks,7) ;peak values
	varv = fltarr(nchunks,7) ;variances

	ylim_start_time=systime(1)
	dmessage="using nchunks in fbk1_7pk: "+string(nchunks)
	dprint,verbose=verbose,dlevel=2,dmessage

	for j=0,6 do begin
		get_data,'rbsp'+probe +'_fbk1_7pk_'+strtrim(j,2),data=goob
		for i=0L,nchunks-1*sz do begin
			maxv[i,j] = max(goob.y[i*sz:i*sz+sz])
			tst = moment(goob.y[i*sz:i*sz+sz])
			varv[i,j] = sqrt(tst[1])
			mom[i,j]=total(goob.y[i*sz:i*sz+sz])/sz
		endfor
	endfor

	stdevv = sqrt(varv)
	ratio = maxv/stdevv


	;Reject values (for the y-scaling) where the max value is >> than the standard deviation
	good = bytarr(nchunks,7)
	good[*] = 1b

	for j=0,6 do begin
		goo = where(ratio[*,j] ge maxrat)
		if goo[0] ne -1 then good[goo,j] = 0b
	endfor

	good = float(good)

	for j=0,6 do ylim,'rbsp'+probe +'_fbk1_7pk_'+strtrim(j,2),0.,max(maxv[*,j]*good[*,j]),0
	if keyword_set(combine) then for j=0,6 do ylim,'rbsp'+probe +'_fbk1_7comb_'+strtrim(j,2),0,max(maxv[*,j]*good[*,j]),0

	dmessage='ylim fbk1_7pk runtime (s): '+string(systime(1)-ylim_start_time)
	dprint,verbose=verbose,dlevel=2,dmessage

	;Keep only 1 number after decimal place. More just clutters up plot
	for j=0,6 do options,'rbsp'+probe +'_fbk1_7pk_'+strtrim(j,2),'ytickformat','(f6.1)'
endif


;---------
;FBK2 - 7
;---------

tst = tnames('rbsp'+probe +'_fbk2_7pk_0')
if tst ne '' then begin
	get_data,'rbsp'+probe +'_fbk2_7pk_0',data=goo0

	nchunks = floor(n_elements(goo0.y)/sz)
	mom = fltarr(nchunks,7)
	maxv = fltarr(nchunks,7) ;peak values
	varv = fltarr(nchunks,7) ;variances

	ylim_start_time=systime(1)
	dmessage="using nchunks in fbk2_7pk: "+string(nchunks)
	dprint,verbose=verbose,dlevel=2,dmessage

	for j=0,6 do begin
		get_data,'rbsp'+probe +'_fbk2_7pk_'+strtrim(j,2),data=goob
		for i=0L,nchunks-1*sz do begin
			maxv[i,j] = max(goob.y[i*sz:i*sz+sz])
			tst = moment(goob.y[i*sz:i*sz+sz])
			varv[i,j] = sqrt(tst[1])
			mom[i,j]=total(goob.y[i*sz:i*sz+sz])/sz
		endfor
	endfor

	stdevv = sqrt(varv)
	ratio = maxv/stdevv


	;Reject values (for the y-scaling) where the max value is >> than the standard deviation
	good = bytarr(nchunks,7)
	good[*] = 1b

	for j=0,6 do begin
		goo = where(ratio[*,j] ge maxrat)
		if goo[0] ne -1 then good[goo,j] = 0b
	endfor

	good = float(good)

	for j=0,6 do ylim,'rbsp'+probe +'_fbk2_7pk_'+strtrim(j,2),0.,max(maxv[*,j]*good[*,j]),0
	if keyword_set(combine) then for j=0,6 do ylim,'rbsp'+probe +'_fbk2_7comb_'+strtrim(j,2),0,max(maxv[*,j]*good[*,j]),0

	dmessage='ylim fbk2_7pk runtime (s): '+string(systime(1)-ylim_start_time)
	dprint,verbose=verbose,dlevel=2,dmessage

	;Keep only 1 number after decimal place. More just clutters up plot
	for j=0,6 do options,'rbsp'+probe +'_fbk2_7pk_'+strtrim(j,2),'ytickformat','(f6.1)'
endif




dmessage='runtime (s): '+string(systime(1)-start_time)
dprint,verbose=verbose,dlevel=2,dmessage

end
