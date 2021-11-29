; rbsp_b1_predict_plot.pro

;Overplot the predicted buffer location over the current buffer location plot
;This is a useful tool for visualizing the future location of the B1 memory pointer location. 
;You can implement B1 collection times, pointer jumps, and indicate regions that you'd
;rather not overwrite. 

;Written by Aaron W Breneman, University of Minnesota, December, 2013


;First need to run:  .run b1_status_crib

;In order to check out exact memory locations at specific times use
    ;ctime, t, buffer_index, npoints=2, /exact
    ;record_rate=(buffer_index[1] - buffer_index[0]) / (t[1] - t[0])
    ;print, 'B1 record rate (blocks/sec):', record_rate

;-----------------------------------------------------------------------------------------

;Reference (memory) location of B1 pointer. All collection times are based on this location. 
;Probably need to update this from time-to-time b/c the prediction times may drift
;from actual collection times. 
;Update using ctime, t, buffer_index, npoints=2, /exact
;Only included predicted times AFTER the reference location

sz = 262144L    ;size of memory buffer (blocks)

cloca = 160415d
clocb = 209006d

;-----------------------------------------------------------------------------------------
; jump id and its start times (AFTER REFERENCE TIME ONLY)
jumpa_s = time_double(['2014-03-08/15:24','2014-03-10/03:16','2014-03-11/15:08'])
nb_jumpa = [248020D,248020D,248020D]

jumpb_s = time_double(['2014-03-08/23:55','2014-03-11/23:46'])
nb_jumpb = [112202D,112202D]

;-----------------------------------------------------------------------------------------

;Start and stop times of requested data playback that we want to be sure to protect (black diamonds are changed to blue diamonds)
tmp = [['2014-02-28/00:35', '2014-02-28/00:35:01']]
tpa0 = time_double(transpose(tmp[0,*]))
tpa1 = time_double(transpose(tmp[1,*]))

tmp =  [['2014-02-22/03:10', '2014-02-22/05:56'],$
        ['2014-02-23/04:00', '2014-02-23/05:00'],$
        ['2014-02-23/07:00', '2014-02-23/08:00']]
tpb0 = time_double(transpose(tmp[0,*]))
tpb1 = time_double(transpose(tmp[1,*]))


;-----------------------------------------------------------------------------------------

;Define date collection rate for each collection time
ratea = [16384,16384,16384,16384,16384,16384,16384,16384,16384,16384,16384,16384,$
    16384,16384,16384,16384,16384,16384,16384,16384]
rateb = [4096,4096,4096,4096,4096,4096,4096,4096,4096,4096,4096,4096,4096,$
    4096,4096,4096,4096,4096,4096,4096,4096]

;Start and end times of collection on A (AFTER REFERENCE TIME ONLY!!!)
tmp = [['2014-03-08/15:25', '2014-03-08/20:25'],$
       ['2014-03-09/00:23', '2014-03-09/05:23'],$
       ['2014-03-09/09:21', '2014-03-09/14:21'],$
       ['2014-03-09/18:19', '2014-03-09/23:19'],$
       ['2014-03-10/03:17', '2014-03-10/08:17'],$
       ['2014-03-10/12:15', '2014-03-10/17:15'],$
       ['2014-03-10/21:13', '2014-03-11/02:13'],$
       ['2014-03-11/06:11', '2014-03-11/11:11'],$
       ['2014-03-11/15:09', '2014-03-11/20:09'],$
       ['2014-03-12/00:07', '2014-03-12/05:07'],$
       ['2014-03-12/09:05', '2014-03-12/14:05'],$
       ['2014-03-12/18:03', '2014-03-12/23:03'],$

       ['2014-03-13/03:01', '2014-03-13/08:01'],$
       ['2014-03-13/11:59', '2014-03-13/16:59'],$
       ['2014-03-13/20:56', '2014-03-14/01:56'],$
       ['2014-03-14/05:54', '2014-03-14/10:54'],$
       ['2014-03-14/14:52', '2014-03-14/19:52'],$
       ['2014-03-14/23:50', '2014-03-15/04:50'],$
       ['2014-03-15/08:48', '2014-03-15/13:48'],$
       ['2014-03-15/17:46', '2014-03-15/22:46']]


timea_s = time_double(transpose(tmp[0,*]))
timea_e = time_double(transpose(tmp[1,*]))

;Start and end times of collection on B (AFTER REFERENCE TIME ONLY!!!)
tmp = [['2014-03-08/05:59', '2014-03-08/11:59'],$
       ['2014-03-08/14:58', '2014-03-08/20:58'],$
       ['2014-03-08/23:56', '2014-03-09/05:56'],$
       ['2014-03-09/08:55', '2014-03-09/14:55'],$
       ['2014-03-09/17:54', '2014-03-09/23:54'],$
       ['2014-03-10/02:53', '2014-03-10/08:53'],$
       ['2014-03-10/11:51', '2014-03-10/17:51'],$
       ['2014-03-10/20:50', '2014-03-11/02:50'],$
       ['2014-03-11/05:49', '2014-03-11/11:49'],$
       ['2014-03-11/14:48', '2014-03-11/20:48'],$
       ['2014-03-11/23:47', '2014-03-12/05:47'],$
       ['2014-03-12/08:45', '2014-03-12/14:45'],$
       ['2014-03-12/17:44', '2014-03-12/23:44'],$

       ['2014-03-13/02:43', '2014-03-13/08:43'],$
       ['2014-03-13/11:42', '2014-03-13/17:42'],$
       ['2014-03-13/20:41', '2014-03-14/02:41'],$
       ['2014-03-14/05:40', '2014-03-14/11:40'],$
       ['2014-03-14/14:38', '2014-03-14/20:38'],$
       ['2014-03-14/23:37', '2014-03-15/05:37'],$
       ['2014-03-15/08:36', '2014-03-15/14:36'],$
       ['2014-03-15/17:35', '2014-03-15/23:35']]

timeb_s = time_double(transpose(tmp[0,*]))
timeb_e = time_double(transpose(tmp[1,*]))


;Do some quick array size checks
if n_elements(ratea) ne n_elements(timea_s) then stop
if n_elements(ratea) ne n_elements(timea_e) then stop
if n_elements(rateb) ne n_elements(timeb_s) then stop
if n_elements(rateb) ne n_elements(timeb_e) then stop
if n_elements(timea_s) ne n_elements(timea_e) then stop
if n_elements(timeb_s) ne n_elements(timeb_e) then stop
if n_elements(jumpa_s) ne n_elements(nb_jumpa) then stop
if n_elements(jumpb_s) ne n_elements(nb_jumpb) then stop
if n_elements(tpa0) ne n_elements(tpa1) then stop
if n_elements(tpb0) ne n_elements(tpb1) then stop


;Blocks per second
ratea2 = ratea/(16000./3.)
rateb2 = rateb/(16000./3.)


;Number of blocks to record for each collection time (3 blocks/sec for 16K)
hopva = (timea_e - timea_s)*ratea2
hopvb = (timeb_e - timeb_s)*rateb2

jumpva = nb_jumpa
jumpvb = nb_jumpb

jump_or_collect_a = [replicate(1,n_elements(timea_s)),replicate(2,n_elements(jumpa_s))]
jump_or_collect_b = [replicate(1,n_elements(timeb_s)),replicate(2,n_elements(jumpb_s))]


;combine collection times and jump times
timea_s = [timea_s,jumpa_s]
timeb_s = [timeb_s,jumpb_s]
timea_e = [timea_e,jumpa_s+0.1]
timeb_e = [timeb_e,jumpb_s+0.1]
jumpatime = [jumpa_s,jumpa_s+0.1]
jumpbtime = [jumpb_s,jumpb_s+0.1]
skipva = [hopva,jumpva]
skipvb = [hopvb,jumpvb]


sta = sort(timea_s)
stb = sort(timeb_s)

joc_a = jump_or_collect_a[sta]
joc_b = jump_or_collect_b[stb]


;Final sorted start and end times for each collection and jump
timea_s = timea_s[sta]
timea_e = timea_e[sta]
timeb_s = timeb_s[stb]
timeb_e = timeb_e[stb]
;Final sorted values of the buffer hop and jump for each collection and jump
incrementva = skipva[sta]
incrementvb = skipvb[stb]

;-------------------------------------------------------------------------------
;Print out timeline of collection and jumps

rr=0
print,'-----------------------------------------------------------------------------------'
print,'RBSP-A'
for uu=0,n_elements(timea_s)-1 do begin

    output = time_string(timea_s[uu])+' to '+ time_string(timea_e[uu])
    if joc_a[uu] eq 2 then output = 'JUMP    ' + output
    if joc_a[uu] eq 1 then output = 'COLLECT ' + output+' rate='+strtrim(floor(ratea[rr]),2)+ ' S/s|'
    if joc_a[uu] eq 1 then output += ' ' + strtrim(floor((timea_e[uu] - timea_s[uu])/60.),2) + ' minutes|'
    if joc_a[uu] eq 1 then output += ' ' + strtrim(floor((timea_e[uu] - timea_s[uu])*ratea2[rr]),2) + ' blocks'
    if joc_a[uu] eq 1 then rr++
    print,output
endfor
print,'-----------------------------------------------------------------------------------'
rr=0
print,'RBSP-B'
for uu=0,n_elements(timeb_s)-1 do begin

    output = time_string(timeb_s[uu])+' to '+ time_string(timeb_e[uu])
    if joc_b[uu] eq 2 then output = 'JUMP    ' + output
    if joc_b[uu] eq 1 then output = 'COLLECT ' + output+' rate='+strtrim(floor(rateb[rr]),2)+ ' S/s|'
    if joc_b[uu] eq 1 then output += ' ' + strtrim(floor((timeb_e[uu] - timeb_s[uu])/60.),2) + ' minutes|'
    if joc_b[uu] eq 1 then output += ' ' + strtrim(floor((timeb_e[uu] - timeb_s[uu])*rateb2[rr]),2) + ' blocks'
    if joc_b[uu] eq 1 then rr++
    print,output
endfor
print,''
;-------------------------------------------------------------------------------

;Future memory locations
mema = dblarr(n_elements(incrementva))
memb = dblarr(n_elements(incrementvb))

; check first location, if it's jump, then do nothing, otherwise use start location.
idx = where(timea_s[0] eq jumpa_s, cnt)
if cnt ne 0 then mema[0] = incrementva[0] else mema[0] = cloca+incrementva[0]
for i=1,n_elements(incrementva)-1 do begin
    idx = where(timea_s[i] eq jumpa_s, cnt)
    ; if it's jump, then do nothing, otherwise increment.
    if cnt ne 0 then mema[i] = incrementva[i] else mema[i] = mema[i-1]+incrementva[i]
endfor

idx = where(timeb_s[0] eq jumpb_s, cnt)
if cnt ne 0 then mema[0] = incrementvb[0] else memb[0] = clocb+incrementvb[0]
for i=1,n_elements(incrementvb)-1 do begin
    idx = where(timeb_s[i] eq jumpb_s, cnt)
    if cnt ne 0 then memb[i] = incrementvb[i] else memb[i] = memb[i-1]+incrementvb[i]
endfor

;Take into account circular nature of buffer
mema = mema mod sz
memb = memb mod sz
mema = floor(mema)
memb = floor(memb)

memas = shift(mema,1)
memas[0] = floor(cloca mod sz)
membs = shift(memb,1)
membs[0] = floor(clocb mod sz)

;combine start and end times of each collection interval
memaf = [memas,mema]
membf = [membs,memb]

;Sort chronologically
timea = [timea_s,timea_e]
timeb = [timeb_s,timeb_e]
sa = sort(timea)
sb = sort(timeb)
timea = timea[sa]
timeb = timeb[sb]
memaf = memaf[sa]
membf = membf[sb]

;Treat the wraparound bug.
memaf = float(memaf)    ; float because need !values.f_nan.
nmemf = n_elements(memaf)
i = 1
while i lt nmemf do begin
    idx = where(timea[i] eq jumpatime, cnt)
    if cnt ne 0 or memaf[i] ge memaf[i-1] then begin
        i+=1
        continue    ; jump, or normal situation, otherwise overflow.
    endif
    wraptime = timea[i-1]+(timea[i]-timea[i-1])*(sz-memaf[i-1])/(memaf[i]+sz-memaf[i-1])
    memaf = [memaf[0:i-1],sz-1,!values.f_nan,0,memaf[i:*]]   ; add f_nan to eliminate vertical line from sz to 0.
    timea = [timea[0:i-1],wraptime-1,wraptime,wraptime+1,timea[i:*]]
    i+=3
    nmemf = n_elements(memaf)
endwhile

membf = float(membf)
nmemf = n_elements(membf)
i = 1
while i lt nmemf do begin
    idx = where(timeb[i] eq jumpbtime, cnt)
    if cnt ne 0 or membf[i] ge membf[i-1] then begin
        i+=1
        continue    ; jump, or normal situation, otherwise overflow.
    endif
    wraptime = timeb[i-1]+(timeb[i]-timeb[i-1])*(sz-membf[i-1])/(membf[i]+sz-membf[i-1])
    membf = [membf[0:i-1],sz-1,!values.f_nan,0,membf[i:*]]
    timeb = [timeb[0:i-1],wraptime-1,wraptime,wraptime+1,timeb[i:*]]
    i+=3
    nmemf = n_elements(membf)
endwhile

;Create a tplot variable with the future memory locations
store_data,'future_a',data={x:timea,y:memaf}
store_data,'future_b',data={x:timeb,y:membf}
options,['future_a','future_b'],'colors',250
options,['future_a','future_b'],'thick',3


;Create a tplot variable with horizontal lines to represent the jumped memory locations


;Find value at jump location
get_data,'future_a',data=bia
get_data,'future_b',data=bib


;use last element of the jump array
gooa = where(bia.x ge jumpa_s[n_elements(jumpa_s)-1])
goob = where(bib.x ge jumpb_s[n_elements(jumpb_s)-1])


v0a = bia.y[gooa[0]]
v1a = v0a + nb_jumpa[n_elements(jumpa_s)-1]
v1a = v1a mod sz

v0b = bib.y[goob[0]]
v1b = v0b + nb_jumpb[n_elements(jumpb_s)-1]
v1b = v1b mod sz

t0 = time_double('2012-01-01/00:00')
t1 = time_double('2050-01-01/00:00')

store_data,'jump_a1',data={x:[t0,t1],y:[v0a,v0a]}
store_data,'jump_a2',data={x:[t0,t1],y:[v1a,v1a]}
store_data,'jump_b1',data={x:[t0,t1],y:[v0b,v0b]}
store_data,'jump_b2',data={x:[t0,t1],y:[v1b,v1b]}

get_data,'rbspa_efw_b1_fmt_block_index2',data=gootmpa
get_data,'rbspb_efw_b1_fmt_block_index2',data=gootmpb
gootmpa2 = gootmpa
gootmpb2 = gootmpb

gootmpa2.y = !values.f_nan
gootmpb2.y = !values.f_nan



for vv=0,n_elements(tpa0)-1 do begin
    boob = where((gootmpa.x ge tpa0[vv]) and (gootmpa.x le tpa1[vv]))
    if boob[0] ne -1 then gootmpa2.y[boob] = gootmpa.y[boob]
endfor
for vv=0,n_elements(tpb0)-1 do begin
    boob = where((gootmpb.x ge tpb0[vv]) and (gootmpb.x le tpb1[vv]))
    if boob[0] ne -1 then gootmpb2.y[boob] = gootmpb.y[boob]
endfor



store_data,'rbspa_efw_b1_fmt_block_index3',data=gootmpa2
store_data,'rbspb_efw_b1_fmt_block_index3',data=gootmpb2
options,'rbsp?_efw_b1_fmt_block_index3','colors',100
options,'rbsp?_efw_b1_fmt_block_index3','psym',4


store_data,'comba',data=['rbspa_efw_b1_fmt_block_index_cutoff','rbspa_efw_b1_fmt_block_index','rbspa_efw_b1_fmt_block_index2','rbspa_efw_b1_fmt_block_index3','future_a','jump_a1','jump_a2']
store_data,'combb',data=['rbspb_efw_b1_fmt_block_index_cutoff','rbspb_efw_b1_fmt_block_index','rbspb_efw_b1_fmt_block_index2','rbspb_efw_b1_fmt_block_index3','future_b','jump_b1','jump_b2']



ylim,['comba','combb'],0,sz

options,'rbsp?_b1_status','panel_size',0.5
tplot,['comba','rbspa_b1_status','combb','rbspb_b1_status']
timebar,jumpa_s,color=50,varname=['comba','rbspa_b1_status']
timebar,jumpb_s,color=50,varname=['combb','rbspb_b1_status']

print,'type .c to print the plot to the desktop'
stop

pcharsize_saved=!p.charsize
pfont_saved=!p.font
pcharthick_saved=!p.charthick
pthick_saved=!p.thick

set_plot,'Z'
rbsp_efw_init,/reset ; try to get decent colors in the Z buffer

device,set_resolution=[3200,2400],set_font='helvetica',/tt_font,set_character_size=[28,35]

!p.thick=4.
!p.charthick=4.

options,['comba','combb'],'ytickformat','(I6.6)'

tplot_options,'xmargin',[14,12]
tplot
timebar,jumpa_s,color=50,varname=['comba','rbspa_b1_status']
timebar,jumpb_s,color=50,varname=['combb','rbspb_b1_status']


; take snapshot of z buffer
snapshot=tvrd()
device,/close

; convert snapshot from index colors to true colors
tvlct,r,g,b,/get

sz=size(snapshot,/dimensions)
snapshot3=bytarr(3,sz[0],sz[1])
snapshot3[0,*,*]=r[snapshot]
snapshot3[1,*,*]=g[snapshot]
snapshot3[2,*,*]=b[snapshot]

; shrink snapshot
xsize=800
ysize=600
snapshot3=rebin(snapshot3,3,xsize,ysize)

print, 'saving png ...'
; write a png
write_png,'~/Desktop/b1_status_predict.png',snapshot3

set_plot,'X'
rbsp_efw_init,/reset
!p.charsize=pcharsize_saved
!p.font=pfont_saved
!p.charthick=pcharthick_saved
!p.thick=pthick_saved


end
