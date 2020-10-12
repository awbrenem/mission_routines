 PRO cxform, oldvec, newvec, xmatrix;;  transforms a vector "oldvec" into a vector "newvec" in another 
;	coordinate system.  xmatrix(na,nc) is a set of vectors for the 
;	axes of the new system, expressed in the old system.  na is
;	obviously the number of the axis;
  version = '15 Feb 2009'
;  Version = '18 Aug 2011'  Version = '25 Aug 2011'; changed 15 Feb so that second index is axis number, to agree with eigenQL
; changed 19 Aug 2011.  Check of indices, see TDS_3Wdecay.pro, search Testmatrix;	previously wrong; 25 Aug 2011, check against Z par B coordinate system, OK; ;
for na = 0,2 do begin
  newvec(na) = xmatrix[0,na]*oldvec[0] + xmatrix[1,na]*oldvec[1] + $
	xmatrix[2,na]*oldvec[2]
endfor
;print,'cx',oldvec[0],newvec[0],xmatrix[*,0]
;print,'cx',oldvec[1],newvec[1],xmatrix[*,1]
;print,'cx',oldvec[2],newvec[2],xmatrix[*,2]
;stopreturnend