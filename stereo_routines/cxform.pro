 PRO cxform, oldvec, newvec, xmatrix



  version = '15 Feb 2009'


for na = 0,2 do begin
  newvec(na) = xmatrix[0,na]*oldvec[0] + xmatrix[1,na]*oldvec[1] + $
	xmatrix[2,na]*oldvec[2]
endfor



