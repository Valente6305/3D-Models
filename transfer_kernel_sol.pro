pro transfer_kernel_sol,ress, Lsource, transf, Dsource=Dsource,  $
              dist=dist, print=print, reset=reset, tau=tau

; J.E.R. Costa 02 Janeiro 2010
; Calculos realizados baseado no formalismo apresentado em
; Ramaty (programa: versao Paulo)
;
; res e' a estrutura resultante do gyro_trans_sol
; res carrega os coeficientes de emissão e auto absorção aord/axtr e eord/extr
; transf e' a estrutura resultante
; Dsource e' o diametro da fonte em arcsec
; LSource = profundidade da fonte em cm
;===================================================================
; nf=n_elements(res.freq)
nn=n_elements(ress.res.h)
res=ress.res
;===================================================================
; Constantes e Equacoes caracteristicas do Plasma:

if not keyword_set(dsource) then dsource=max(Lsource); as vezes L e'vect
if not keyword_set(reset) then reset=0
tth = dsource*!dtor/3600.  ; Theta em radianos.(dsource must be arcsec!!)
AU = 1.49597870d13         ; Unidade Astronômica em cm
if keyword_set(dist) then distance=dist else distance=AU
cgs2sfu = 1.d19

; Tamanho da fonte: [cilindro -> voxel]
Vol = LSource * (tth * distance)^2d         ; volume de um cubo
arc2cm = 1d0 * !dtor/3600d0 * 1.49597870d13 ; arcsec para cm no sol
m0 = 9.1094d-28                             ; electron mass [g]
c = 2.998d10                                ; [cm/s] speed of light
e = 4.803d-10                               ; electron charge
kB= 1.3806d-16                              ; Boltzmann cte (erg /deg)
E0 = m0 * c^2.                              ; Electron rest mass em erg
E0k = E0 / 1.6022d-12 / 1d3                 ; em keV

wp=double(res.w_p)
wc=double(res.w_gyro)

if (size(LSource))[0] eq 0 then s=double(LSource) else stop
                        
if total((size(transf))[0:1]) eq 0 or reset then begin
   Io=replicate(0.d0,nn)
   Ix=replicate(0.d0,nn)
endif else begin
   Io=transf.Io ;[0]
   Ix=transf.Ix ;[0]
endelse

tau=fltarr(nn,2)
;=========================================================
; Construindo a transferencia parcial ate' o s do momento


angsol = tth^2.
tau[*,0]=res.aord*s
tau[*,1]=res.axtr*s
pp=where(finite(tau[*,0]) eq 0)
if pp[0] ne -1 then tau[pp,0]=0.
pp=where(finite(tau[*,1]) eq 0)
if pp[0] ne -1 then tau[pp,1]=0.
raz_ord=res.eord/res.aord
raz_ext=res.extr/res.axtr
pp=where(finite(raz_ord) eq 0)
if pp[0] ne -1 then raz_ord[pp]=0.
pp=where(finite(raz_ext) eq 0)
if pp[0] ne -1 then raz_ext[pp]=0.

;Iot=raz_ord*(1.0-exp(-tau[*,0]))+Io*exp(-tau[*,0]) ; comentei essas duas
;Ixt=raz_ext*(1.0-exp(-tau[*,1]))+Ix*exp(-tau[*,1])

;====================&&&&&&&&&&&===================
;; RADIATIVE TRANSFER
 Iot=dblarr(nn)
 Ixt=dblarr(nn)
 tau_o=res.aord*s
 tau_x=res.axtr*s

 toLow = where(tau_o LT 0.001)
 toMed = where((tau_o GE 0.001) AND (tau_o LE 20.))
 toHig = where(tau_o GT 20.)
 txLow = where(tau_x LT 0.001)
 txMed = where((tau_x GE 0.001) AND (tau_x LE 20.))
 txHig = where(tau_x GT 20.)
;stop
 IF toLow[0] NE -1 THEN Iot[toLow] = res[toLow].eord * s
 IF toMed[0] NE -1 THEN Iot[toMed] = raz_ord[toMed]*(1d0-exp(-tau_o[toMed]))
 IF toHig[0] NE -1 THEN Iot[toHig] = raz_ord[toHig] 
 IF txLow[0] NE -1 THEN Ixt[txLow] = res[txLow].extr * s
 IF txMed[0] NE -1 THEN Ixt[txMed] = raz_ext[txMed]*(1d0-exp(-tau_x[txMed]))
 IF txHig[0] NE -1 THEN Ixt[txHig] = raz_ext[txHig]

Iot+=Io*exp(-tau_o) 
Ixt+=Ix*exp(-tau_x)


;====================************==================
;Iot=Iot*angsol
;Ixt=Ixt*angsol
;==========================================================

trt={w:0d,h:0d,angle:0d,Io:0d,Ix:0d,$
       o:0d,cgs2sfu:0d, $
       L_total:0d,B:0d,Nel:0d,Omega:0d}
transf=replicate(trt,nn)
transf.h=res.h                         ; n# harmonico
transf.angle=res.angle                 ; ang. de visada
transf.Io=Iot                          ; emissiv. ordinaria
transf.Ix=Ixt                          ; emissiv. extraordinaria
transf.o=tth^2                         ; ang. solido em seg^2
transf.cgs2sfu=replicate(cgs2sfu,nn)   ; fator de conversao p/ sfu
transf.L_total=Lsource                 ; profundidade da fonte
transf.B=res.b                         ; Inducao magn.
transf.Nel=res.nel                     ; dens. n# elet. nao termicos
transf.Omega=tth                       ; diametro angular

end
