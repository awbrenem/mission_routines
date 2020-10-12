;Station locations at: http://wdc.kugi.kyoto-u.ac.jp/aedir/ae2/AEObs.html

;Get data from individual stations that comprise the AE network from:
;http://www.intermagnet.org/data-donnee/dataplot-eng.php?

;TABLE 1 - List of AE(12) Stations.
;IAGA	   Geographic Coord.	  Geomagnetic Coord.
;Observatory	Code	       Lat.(°N)	Long.(°E)	Lat.(°N)	Long.(°E)
;Abisko	ABK	             68.36	  18.82	    66.04	    115.08
;Dixon Island	DIK	       73.55	  80.57	    63.02	    161.57
;Cape Chelyuskin	CCS	     77.72	  104.28	  66.26	    176.46
;Tixie Bay	TIK	           71.58	  129.00	  60.44	    191.41
;Cape Wellen	CWE	         66.17	  190.17	  61.79	    237.10
;Barrow	BRW	             71.30	  203.25	  68.54	    241.15
;College	CMO	             64.87	  212.17	  64.63	    256.52
;Yellowknife	YKC	         62.40    245.60	  69.00	    292.80
;Fort Churchill	FCC	     58.80	  265.90	  68.70	    322.77
;Poste-de-la-Baleine	PBQ  55.27	  282.22	  66.58	    347.36
;Narsarsuaq NAQ	         61.20	  314.16	  71.21	    36.79
;Leirvogur	LRV	           64.18	  338.30	  70.22	    71.04


pro read_ae_station_data


path = '/Users/aaronbreneman/Desktop/code/Aaron/github.umn.edu/AE_magnetometer_plots/data20190425183342/'
fn = 'aaa20140110vmin.min'

;DATE       TIME         DOY     ABKX      ABKY      ABKZ      ABKF   |
;2014-01-10 00:00:00.000 010     11350.90   1569.20  51854.30  53105.70

openr,lun,path+fn,/get_lun
jnk = ''
for i=0,20 do readf,lun,jnk






close,lun





end
