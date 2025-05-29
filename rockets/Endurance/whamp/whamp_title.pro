function whamp_title, plasma, wna
beta = plasma_beta(plasma) & beta = beta[0:1]
title='WNA = '+number_to_string(wna,format='(f12.2)') + '!uo!n'

lab1 = ' '+ strmid(plasma.label,0, strpos(plasma.label,'nT') + 2)
lab2 = ' PLASMA ' + symbol0('beta') +'=' + str_flatten(number_to_string(beta))
ii=where(plasma.species.frac ne 0)
lab3 = plasma.labels[ii[0]] + '!c'
for i=1,n_elements(ii)-2 do lab3 = lab3 + plasma.labels[ii[i]] + '!c'
lab3=lab3 + plasma.labels[ii[i]]
title= title + lab1 + lab2 + '!c' + lab3
return, title
end
