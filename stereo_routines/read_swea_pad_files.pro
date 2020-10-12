;Read STEREO SWEA PAD files for a single day. Output a tplot save file with all
;relevant tplot variables.

;These are the formal calibrated, vetted files the SWEA team wants you to use.
;Find at https://stereo-ssc.nascom.nasa.gov/data/ins_data/impact/level1/ahead/swea/ascii/PAD/
;Unfortunately they're not available in CDF.

;They do have the full distributions in CDF, but these are not properly calibrated
;and include the bad <50 eV data, etc..


path = '/Users/aaronbreneman/Desktop/Research/OTHER/Stuff_for_other_people/Cattell_Cindy/2019_paper_plots/'
;fn = 'STA_L2_SWEA_PAD_20170324_V04.cef.gz'
;fn = 'STA_L2_SWEA_PAD_20170325_V04.cef.gz'
fn = 'STA_L2_SWEA_PAD_20170326_V04.cef.gz'

pre = 'sta_'
filename = pre + '20170326pads'

openr,lun,path+fn,/get_lun
jnk = ''
while jnk ne 'DATA_UNTIL = EOF' do readf,lun,jnk

;Pitch angles
patable = [7.50, 22.50, 37.50, 52.50, 67.50, 82.50, 97.50, 112.50, 127.50, 142.50, 157.50, 172.50]


times = ''
qi0 = ''
qi1 = ''
duration = ''
energy_table = replicate('',16)
delta_plus_energy_table = replicate('',16)
delta_minus_energy_table = replicate('',16)
data2d1 = replicate('',16) & data2d2 = replicate('',16) & data2d3 = replicate('',16)
data2d4 = replicate('',16) & data2d5 = replicate('',16) & data2d6 = replicate('',16)
data2d7 = replicate('',16) & data2d8 = replicate('',16) & data2d9 = replicate('',16)
data2d10 = replicate('',16) & data2d11 = replicate('',16) & data2d12 = replicate('',16)


while not eof(lun) do begin

  dat = ''
  readf,lun,dat
  datv = strsplit(dat,',',/extract)

  times = [times,datv[0]]
  qi0 = [qi0,datv[1]]
  qi1 = [qi1,datv[2]]
  duration = [duration,datv[3]]
  energy_table = [[energy_table],[datv[4:19]]]
  delta_plus_energy_table = [[delta_plus_energy_table],[datv[20:35]]]
  delta_minus_energy_table = [[delta_minus_energy_table],[datv[36:51]]]

  tmp = datv[(16*0)+52:(16*0)+67] & data2d1 = [[data2d1],[tmp]]
  tmp = datv[(16*1)+52:(16*1)+67] & data2d2 = [[data2d2],[tmp]]
  tmp = datv[(16*2)+52:(16*2)+67] & data2d3 = [[data2d3],[tmp]]
  tmp = datv[(16*3)+52:(16*3)+67] & data2d4 = [[data2d4],[tmp]]
  tmp = datv[(16*4)+52:(16*4)+67] & data2d5 = [[data2d5],[tmp]]
  tmp = datv[(16*5)+52:(16*5)+67] & data2d6 = [[data2d6],[tmp]]
  tmp = datv[(16*6)+52:(16*6)+67] & data2d7 = [[data2d7],[tmp]]
  tmp = datv[(16*7)+52:(16*7)+67] & data2d8 = [[data2d8],[tmp]]
  tmp = datv[(16*8)+52:(16*8)+67] & data2d9 = [[data2d9],[tmp]]
  tmp = datv[(16*9)+52:(16*9)+67] & data2d10 = [[data2d10],[tmp]]
  tmp = datv[(16*10)+52:(16*10)+67] & data2d11 = [[data2d11],[tmp]]
  tmp = datv[(16*11)+52:(16*11)+67] & data2d12 = [[data2d12],[tmp]]

endwhile
close,lun & free_lun,lun


;Remove leading element
n = n_elements(times)-1

times = time_double(times[1:n])
qi0 = byte(reform(qi0[1:n]))
qi1 = byte(reform(qi1[1:n]))
duration = float(duration[1:n])
energy_table = float(energy_table[*,1:n])
delta_plus_energy_table = float(delta_plus_energy_table[*,1:n])
delta_minus_energy_table = float(delta_minus_energy_table[*,1:n])
data2d1 = float(transpose(data2d1[*,1:n])) & data2d2 = float(transpose(data2d2[*,1:n]))
data2d3 = float(transpose(data2d3[*,1:n])) & data2d4 = float(transpose(data2d4[*,1:n]))
data2d5 = float(transpose(data2d5[*,1:n])) & data2d6 = float(transpose(data2d6[*,1:n]))
data2d7 = float(transpose(data2d7[*,1:n])) & data2d8 = float(transpose(data2d8[*,1:n]))
data2d9 = float(transpose(data2d9[*,1:n])) & data2d10 = float(transpose(data2d10[*,1:n]))
data2d11 = float(transpose(data2d11[*,1:n])) & data2d12 = float(transpose(data2d12[*,1:n]))


pad_data = [[[data2d1]],[[data2d2]],[[data2d3]],[[data2d4]],[[data2d5]],$
[[data2d6]],[[data2d7]],[[data2d8]],[[data2d9]],[[data2d10]],[[data2d11]],[[data2d12]]]

energies = reform(energy_table[*,0])
store_data,pre+'energyvals',times,transpose(energy_table)
store_data,pre+'pad',data={x:times,y:pad_data,v1:energies,v2:patable}


; set up colors
device,decomposed=0
loadct,39
options,'*pad_*','no_interp',1  ; turn off interpolation




;Pitch angles
;      7.50000      22.5000      37.5000      52.5000      67.5000      82.5000
;      97.5000      112.500      127.500      142.500      157.500      172.500
store_data,pre+'pad_7.5deg',data={x:times,y:reform(pad_data[*,0:7,0]),v:energies[0:7]}
store_data,pre+'pad_22.5deg',data={x:times,y:reform(pad_data[*,0:7,1]),v:energies[0:7]}
store_data,pre+'pad_37.5deg',data={x:times,y:reform(pad_data[*,0:7,2]),v:energies[0:7]}
store_data,pre+'pad_52.5deg',data={x:times,y:reform(pad_data[*,0:7,3]),v:energies[0:7]}
store_data,pre+'pad_67.5deg',data={x:times,y:reform(pad_data[*,0:7,4]),v:energies[0:7]}
store_data,pre+'pad_82.5deg',data={x:times,y:reform(pad_data[*,0:7,5]),v:energies[0:7]}
store_data,pre+'pad_97.5deg',data={x:times,y:reform(pad_data[*,0:7,6]),v:energies[0:7]}
store_data,pre+'pad_112.5deg',data={x:times,y:reform(pad_data[*,0:7,7]),v:energies[0:7]}
store_data,pre+'pad_127.5deg',data={x:times,y:reform(pad_data[*,0:7,8]),v:energies[0:7]}
store_data,pre+'pad_142.5deg',data={x:times,y:reform(pad_data[*,0:7,9]),v:energies[0:7]}
store_data,pre+'pad_157.5deg',data={x:times,y:reform(pad_data[*,0:7,10]),v:energies[0:7]}
store_data,pre+'pad_172.5deg',data={x:times,y:reform(pad_data[*,0:7,11]),v:energies[0:7]}


options,pre+'pad_*deg','spec',0
ylim,pre+'pad_*deg',0.001,10000.,1
tplot,pre+'pad_*deg'


;Energies
;1716.84 1056.95 650.70 400.60 246.62 151.83 93.47 57.54 35.43 21.81 13.43 8.27 5.09 3.13 1.93 1.19
store_data,pre+'pad_1720eV',data={x:times,y:reform(pad_data[*,0,*]),v:patable}
store_data,pre+'pad_1056eV',data={x:times,y:reform(pad_data[*,1,*]),v:patable}
store_data,pre+'pad_650eV',data={x:times,y:reform(pad_data[*,2,*]),v:patable}
store_data,pre+'pad_400eV',data={x:times,y:reform(pad_data[*,3,*]),v:patable}
store_data,pre+'pad_246eV',data={x:times,y:reform(pad_data[*,4,*]),v:patable}
store_data,pre+'pad_152eV',data={x:times,y:reform(pad_data[*,5,*]),v:patable}
store_data,pre+'pad_93eV',data={x:times,y:reform(pad_data[*,6,*]),v:patable}
store_data,pre+'pad_58eV',data={x:times,y:reform(pad_data[*,7,*]),v:patable}


options,pre+'pad_*eV','spec',1
ylim,pre+'pad_*eV',0,180
zlim,pre+'pad_1720eV',1d-5,1d-2,1
zlim,pre+'pad_1056eV',1d-4,1d-1,1
zlim,pre+'pad_650eV',1d-3,1d2,1
zlim,pre+'pad_400eV',1d-2,1d2,1
zlim,pre+'pad_246eV',1d-1,1d2,1
zlim,pre+'pad_152eV',1d0,1d3,1
zlim,pre+'pad_93eV',1d2,1d4,1
zlim,pre+'pad_58eV',1d3,1d4,1
tplot,pre+'pad_*eV'


savevars = [tnames(pre+'pad*')]
tplot_save,savevars,filename='~/Desktop/'+filename


stop

end
