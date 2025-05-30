;+
; Read RBSP keV proton flux.
; Save as rbspx_kev_ion_flux.
;
; set pitch_angle to load data for a specific pitch angle, otherwise load all pitch angles.
;-
pro rbsp_read_kev_proton, time, probe=probe, errmsg=errmsg, pitch_angle=pitch_angle, energy=energy

    pre0 = 'rbsp'+probe+'_'

    ; read 'rbspx_kev_h_flux'
    rbsp_read_mageis, time, id='l3%ion', probe=probe, errmsg=errmsg
    if errmsg ne 0 then return

    var = pre0+'kev_h_flux'
    enbins = get_var_data('FPDU_Energy')
    nenbin = n_elements(enbins)
    enidx = where(finite(enbins), nenbin)
    enbins = enbins[enidx]

    get_data, 'FPDU', uts, dat
    dat = reform(dat[*,enidx,*])>1

    ; apply energy range.
    if n_elements(energy) eq 0 then enidx = findgen(nenbin) else begin
        case n_elements(energy) of
            1: begin
                enidx = where(enbins eq energy, cnt)
                if cnt eq 0 then tmp = min(enbins-energy[0], /absolute, enidx)
            end
            2: begin
                enidx = where(enbins ge energy[0] and enbins le energy[1], cnt)
                if cnt eq 0 then begin
                    errmsg = 'no energy in given range ...'
                    return
                endif
            end
            else: begin
                errmsg = 'wrong # of energy info ...'
                return
            end
        endcase
    endelse
    dat = dat[*,enidx,*]
    enbins = enbins[enidx]
    nenbin = n_elements(enbins)


    ; filter pitch angle.
    pabins = get_var_data('FPDU_Alpha')
    npabin = n_elements(pabins)
    if n_elements(pitch_angle) eq 0 then paidx = findgen(npabin) else begin
        case n_elements(pitch_angle) of
            1: begin
                paidx = where(pabins eq pitch_angle, cnt)
                if cnt eq 0 then tmp = min(pabins-pitch_angle[0], /absolute, paidx)
            end
            2: begin
                paidx = where(pabins ge pitch_angle[0] and pabins le pitch_angle[1], cnt)
                if cnt eq 0 then begin
                    errmsg = 'no pitch angle in given range ...'
                    return
                endif
            end
            else: begin
                errmsg = 'wrong # of pitch angle info ...'
                return
            end
        endcase
    endelse
    dat = reform(dat[*,*,paidx])
    pabins = pabins[paidx]
    npabin = n_elements(pabins)

    ; Average pitch angle if no pitch angle info is provided.
    if n_elements(pitch_angle) eq 0 then begin
        dat = total(dat,3,/nan)/npabin
        npabin = 0
    endif

    ; save data.
    if nenbin eq 1 and npabin eq 1 then begin
        store_data, var, uts, dat
        add_setting, var, /smart, {$
            display_type: 'scalar', $
            ylog: 1, $
            unit: '#/cm!U2!N-s-sr-keV', $
            short_name: 'H!U+!N flux '+sgnum2str(sround(pabins))+'deg, '+sgnum2str(sround(enbins))+'keV'}
    endif else if nenbin eq 1 then begin    ; flux vs pitch angle at certain energy.
        store_data, var, uts, dat, pabins
        add_setting, var, /smart, {$
            display_type: 'list', $
            ylog: 1, $
            unit: '#/cm!U2!N-s-sr-keV', $
            value_unit: 'deg', $
            short_name: 'H!U+!N flux '+sgnum2str(sround(enbins))+' keV'}
    endif else if npabin eq 1 then begin    ; flux vs energy at certain pitch angle.
        yrange = 10d^ceil(alog10(minmax(dat)))>1
        store_data, var, uts, dat, enbins
        add_setting, var, /smart, {$
            display_type: 'list', $
            ylog: 1, $
            yrange: yrange, $
            color_table: 52, $
            unit: '#/cm!U2!N-s-sr-keV', $
            value_unit: 'keV', $
            short_name: 'H!U+!N flux '+sgnum2str(sround(pabins))+' deg'}
    endif else begin
        store_data, var, uts, data
        add_setting, var, /smart, {$
            display_type: 'data', $
            ylog: 1, $
            energy: enbins, $
            pitch_angle: pabins}
    endelse


    dt = 10.848
    uniform_time, var, dt

    store_data, 'FPDU'+['','_Alpha','_Energy'], /delete

end

time = time_double(['2014-08-28/03:00','2014-08-28/15:00'])
probe = 'b'
pitch_angle = 90
energy = [100,300]

time = time_double(['2014-08-28/03:00','2014-08-28/15:00'])
probe = 'a'
pitch_angle = 45
energy = [30,300]

rbsp_read_kev_proton, time, probe=probe, pitch_angle=pitch_angle, energy=energy
end
