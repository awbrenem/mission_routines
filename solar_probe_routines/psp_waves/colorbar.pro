; $Id: //depot/Release/ENVI53_IDL85/idl/idldir/lib/graphics/colorbar.pro#1 $
;
; Copyright (c) 2000-2015, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;

;----------------------------------------------------------------------------
;+
; :Description:
;    Create IDL Colourbar graphic.
;
; :Params:
;    
; :Keywords:
;    _REF_EXTRA
;
; :Returns:
;    Object Reference
;-
function colorbar, DEBUG=debug, _REF_EXTRA=_extra
  compile_opt idl2, hidden
@graphic_error

  return, OBJ_NEW('Colorbar', _EXTRA=_extra)

end
