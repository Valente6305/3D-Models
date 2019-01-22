Function compr,x,y, px=px, mm=mm
; ==========================================================
; Esta funcao mede o comprimento de arco (das linhas) desde 
; o primeiro ponto do vetor enviado.
; Se nao houver keywords sai em pixels (original)
; Se enviar px sai em cm
; Se enviar px e mm (max) sai normalizado para o max (mm)
; ==========================================================

AU = 1.49597870d13                   ; Astronomic Unit
arc2cm = !dtor/3600d0 * AU           ; arcsec to cm in Sun
if keyword_set(px) then fact=px*arc2cm else fact=1.0
if keyword_set(mm) then norma=mm*fact else norma=1.0

n=n_elements(x)
compr=fltarr(n)
for i=1,n-2l do begin
   cir_3pnt,x[i-1l:i+1l],y[i-1l:i+1l],r,x0,y0
   dh=r 
   if dh ne 0 then begin 
      dop=sqrt((x[i-1l]-x[i])^2+(y[i-1l]-y[i])^2)/2.0
      compr[i]=2.0*asin(dop/dh)*r  
   endif else compr[i]=sqrt((x[i-1l]-x[i])^2+(y[i-1l]-y[i])^2)
endfor 
; Ultimo ponto (n-1)
if dh ne 0 then begin
    dop=sqrt((x[i-1l]-x[i])^2+(y[i-1l]-y[i])^2)/2.0
    compr[i]=2.0*asin(dop/dh)*r 
 endif else compr[i]=sqrt((x[i-1l]-x[i])^2+(y[i-1l]-y[i])^2)
return,total(compr*fact/norma,/cum)
end
