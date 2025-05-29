
function whamp_fields,q

a = whamp_field(q[0].sol)
for i=1,n_elements(q)-1 do a = [a, whamp_field(q[i].sol)]

return,a
end

