;+
; NAME:cor
;
;
;
; PURPOSE:Setup color for device
;
;
;
; CATEGORY:main programm
;
;
;
; CALLING SEQUENCE:cor
;
;
;
; INPUTS:none
;
;
;
; OPTIONAL INPUTS:none
;
;
;
; KEYWORD PARAMETERS:
;
;
;
; OUTPUTS:red temperature
;
;
;
; OPTIONAL OUTPUTS: loadct,3 {from 1 to 40}
;
;
;
; COMMON BLOCKS:
;
;
;
; SIDE EFFECTS:
;
;
;
; RESTRICTIONS:
;
;
;
; PROCEDURE:
;
;
;
; EXAMPLE:cor
;
;
;
; MODIFICATION HISTORY: Valente Cuambe August,25th 2014
;
;-
pro cor

device,true=24,decompose=0,retain=2
loadct,3

end
