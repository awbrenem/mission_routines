FUNCTION endurance_lfdsp_construct_15bitword, data
;
; Combine two words to full 15-bit
;TM word 29             30              31
;----------------|-------------|----------------
;HHHHHHHHHHbbbbbaaaaaHHHHHHHHHH
;Where W29+aaaaa makes up first channel and W31+bbbbb makes up second channel
;
	print, '     dynamo2_lfdsp_construct_15bitword...'

	shift_word='4000'XL  ; 16,384 dec = 2^14

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; First channel
    ;
    ; high 10 bits at word index 0, 3, 6
	msb=ishft(data[0:*:3,*] AND 1023, 5)
	; low 5 bits at word index 1, 4, 7
	lsb=data[1:*:3,*] AND 31
	
	; Convert arrays to vectors
	; help, msb
	; help, lsb
    nvals_lsb = size(lsb,/n_elements)
    lsb_1d = reform(lsb, 1, nvals_lsb)
    lsb_vec= reform(lsb_1d)
    nvals_msb = size(msb,/n_elements)
    msb_1d = reform(msb, 1, nvals_msb)
    msb_vec= reform(msb_1d)
    ; help, lsb_vec
    ; help, msb_vec
 
	; Add low bits and high bits, then the infamous XOR correction
	word15bit_set1=(lsb_vec+msb_vec) XOR shift_word	
	
    help, word15bit_set1
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
    ; Second channel
    ;
    ; high 10 bits at word index 2, 5, 8
	msb=ishft(data[2:*:3,*] AND 1023, 5)
	; low 5 bits at word index 1, 4, 6
	lsb=ishft(data[1:*:3,*] AND 992, -5)
	
	; Convert arrays to vectors
	; help, msb
	; help, lsb
    nvals_lsb = size(lsb,/n_elements)
    lsb_1d = reform(lsb, 1, nvals_lsb)
    lsb_vec= reform(lsb_1d)
    nvals_msb = size(msb,/n_elements)
    msb_1d = reform(msb, 1, nvals_msb)
    msb_vec= reform(msb_1d)
    ; help, lsb_vec
    ; help, msb_vec
    
	; Add low bits and high bits, then the infamous XOR correction
	word15bit_set2=(lsb_vec+msb_vec) XOR shift_word	
	
	help, word15bit_set2
    ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

    nvals_set1=size(word15bit_set1,/n_elements)
    nvals_set2=size(word15bit_set2,/n_elements)

	if (nvals_set1 NE nvals_set2) then begin
		print, 'Halting, two sets do not match in size'
		stop
	endif
    
    words15bit=ulonarr(2, nvals_set1)
    
    words15bit[0,*]=word15bit_set1
    words15bit[1,*]=word15bit_set2

return, words15bit
end

;+
; NAME:
;   dynamo2_main_lfdsp_vlf_decode.pro
;
;-
pro endurance_tm1_lfdsp_vlf_decode

print, 'endurance_tm1_lfdsp_vlf_decode.pro...'

fileout='47001_TM1_LFDSP_S5_VLF_mvm.sav'
filein=strarr(2)
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Extract 
print,'Extracting VLF12D and VLF34D...'

bl1234=3.212

filein[0]='47001_TM1_LFDSP_S5VLF1234_VLF12D34D.sav'
restore, filein[0]

print, fileout, format='("fileout= ", a-60)'

dv15bit=endurance_lfdsp_construct_15bitword(dvals)

nvals=size(dv15bit[0,*],/n_elements)

dvlf12raw=transpose(float(dv15bit[1,*]))
dvlf34raw=transpose(float(dv15bit[0,*]))

help, dvlf12raw
help, dvlf34raw

; Normalize to -1 to 1
dvlf12raw_normal=float(dvlf12raw)/(16384.0) - 1.0
dvlf34raw_normal=float(dvlf34raw)/(16384.0) - 1.0

print, min(dvlf12raw_normal), max(dvlf12raw_normal), format='("dvlf12raw_normal:", 2(e12.5,1x))'

; Calibrate using cal from Paulo (spreadsheet dated 2-11-2022)
A12=  1.265e-1
B12=  5.5604e-5
dvlf12_volts=A12*dvlf12raw_normal+B12
dvlf12_mvm=1000.0*dvlf12_volts/bl1234

A34=  1.245e-1
B34=  5.1386e-5
dvlf34_volts=A34*dvlf34raw_normal+B34
dvlf34_mvm=1000.0*dvlf34_volts/bl1234

; Time tags, 1 ms
tvlf=dblarr(size(dvlf12_mvm, /n_elements))
dtword=1.0d0/30.0d3
dtframe=1.0d-3
tvlf[0:*:3]=tvals-dtframe
tvlf[1:*:3]=tvals-dtframe+dtword
tvlf[2:*:3]=tvals-dtframe+dtword*2.0d0
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Extract 
print,'Extracting VLF24D and VLF32D...'

bl2432=2.2712

filein[1]='47001_TM1_LFDSP_S5VLF2432_VLF24D32D.sav'
restore, filein[1]

dv15bit=endurance_lfdsp_construct_15bitword(dvals)

nvals=size(dv15bit[0,*],/n_elements)

dvlf24raw=transpose(float(dv15bit[0,*]))
dvlf32raw=transpose(float(dv15bit[1,*]))

help, dvlf24raw
help, dvlf32raw

; Normalize to -1 to 1
dvlf24raw_normal=float(dvlf24raw)/(16384.0) - 1.0
dvlf32raw_normal=float(dvlf32raw)/(16384.0) - 1.0

print, min(dvlf24raw_normal), max(dvlf24raw_normal), format='("dvlf24raw_normal:", 2(e12.5,1x))'

; Calibrate using cal from Paulo (spreadsheet dated 2-11-2022)
A24=  1.234e-1
B24=  5.7427e-5
dvlf24_volts=A24*dvlf24raw_normal+B24
dvlf24_mvm=1000.0*dvlf24_volts/bl2432

A32=  1.249e-1
B32=  5.4683e-5
dvlf32_volts=A32*dvlf32raw_normal+B32
dvlf32_mvm=1000.0*dvlf32_volts/bl2432

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Extract 
print,'Extracting VLF41D and VLF13D...'

bl4113=2.2712

filein[1]='47001_TM1_LFDSP_S5VLF4113_VLF41D13D.sav'
restore, filein[1]

dv15bit=endurance_lfdsp_construct_15bitword(dvals)

nvals=size(dv15bit[0,*],/n_elements)

dvlf41raw=transpose(float(dv15bit[0,*]))
dvlf13raw=transpose(float(dv15bit[1,*]))

help, dvlf41raw
help, dvlf13raw

; Normalize to -1 to 1
dvlf41raw_normal=float(dvlf41raw)/(16384.0) - 1.0
dvlf13raw_normal=float(dvlf13raw)/(16384.0) - 1.0

print, min(dvlf41raw_normal), max(dvlf41raw_normal), format='("dvlf41raw_normal:", 2(e12.5,1x))'

; Calibrate using cal from Paulo (spreadsheet dated 2-11-2022)
A41=  1.253e-1
B41=  -5.279e-4
dvlf41_volts=A41*dvlf41raw_normal+B41
dvlf41_mvm=1000.0*dvlf41_volts/bl4113

A13=  1.266e-1
B13=  4.4484e-5
dvlf13_volts=A13*dvlf13raw_normal+B13
dvlf13_mvm=1000.0*dvlf13_volts/bl4113


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Save	
dataunits='dvlf{12|34|24|32|41|13}_mvm = mV/m, tvlf = Seconds since T-0'
calnote_vlf1234='Counts -> Volts: A*dvlf{12|34}_normal + B, A12=  1.265e-1, B12=  5.5604e-5, A34=1.245e-1, B34=5.1386e-5, boomlength=3.212'
calnote_vlf2432='Counts -> Volts: A*dvlf{24|32}_normal + B, A24=  1.234e-1, B24=  5.7427e-5, A32=1.249e-1, B32=5.4683e-5, boomlength=2.2712'
calnote_vlf4113='Counts -> Volts: A*dvlf{41|13}_normal + B, A41=  1.253e-1, B41=  -5.279e-4, A13=1.266e-1, B13=4.4484e-5, boomlength=2.2712'
samplerate=30000
 
save, tvlf, dvlf12_mvm, dvlf34_mvm, dvlf24_mvm, dvlf32_mvm, dvlf41_mvm, dvlf13_mvm, $
    author, calnote_vlf1234, calnote_vlf2432, calnote_vlf4113, dataunits, flight, format, filein, link, samplerate, $
    t0, timetagmethod, timeunits,filename=fileout
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

end
