;Check for NaN values in the ev struct for psp_waves_save
PRO psp_waves_ev_check,ev
    
    ;wherever there are NaN values, or infinite values, replace these by -1 in the output
    ;file in order to avoid formatting errors
    ;ev.freq[where(~finite(ev.freq[*]))] = -1
    ;ev.bw[where(~finite(ev.bw[*]))] = -1
    ev.df_f[where(~finite(ev.df_f[*]))] = -1
   
END ;procedure psp_waves_ev_check