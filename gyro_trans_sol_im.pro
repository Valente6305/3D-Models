pro gyro_trans_sol_im, Delta, Energy, Nel_v, M, Bmag_v, Np_v, ANGLE_v, Freq, res, $
                   bremss=bremss, px=px, problem=problem, temperature=temperature

if not keyword_set(px) then px=4.785e7 ; cm/1 pixel
if not keyword_set(bremss) then bremss=1 ; com Bremsstrhalung
;if (n_elements(where(Nel_v eq Nel_v[0])) ne n_elements(Nel_v)) then begin
;   print, 'Densidade nao termica tem que ser contante!'
;   stop
;endif
dens_nt=Nel_v[0] 
nn=n_elements(Bmag_v)
rest={B:0d,angle:0d,Nel:0d,h:0d,w_gyro:0d, $
      w_p:0d,eord:0d,extr:0d,aord:0d,axtr:0d}  ;,flux:0d}
ress=replicate(rest,nn)
;res=replicate(resf,n_freq)

   eord=dblarr(nn)
   extr=dblarr(nn)
   aord=dblarr(nn)
   axtr=dblarr(nn)
   fluxo=dblarr(nn,/nozero)
   fluxo[*]=1.

; alterando para freq vetor
; nf=n_elements(freq)

   Bmag_v = Bmag_v*1.0d0        ; Gauss
   f_gyro=2.8e6*Bmag_v          ; em Hz
   freq=freq[0]*1.d0
   LSource=px ; 1.386E8              ; em cm  (O QUE E'ISSO!!!)
   in_gyro={Delta:delta,Energy:energy,M:M,Np:Np_v,bremss:bremss,px:px}

; Constantes e Equacoes caracteristicas do Plasma:
   arc2cm = 1d0 * !dtor/3600d0 * 1.49597870d13 ; arcsec para cm no sol
   m0 = 9.1094d-28                             ; electron mass [g]
   AU = 1.495978e13                            ; Unid. Astron. cm
   c = 2.998d10                                ; [cm/s] speed of light
   e = 4.803d-10                               ; electron charge
   kB= 1.3806d-16                              ; Boltzmann cte (erg /deg)
   E0 = m0 * c^2.                              ; Electron rest mass em erg
   E0k = E0 / 1.6022d-12 / 1d3                 ; em keV
   w_gyro = 2.0*!dpi*f_gyro                    ; w_gyro = e * B / (m0 * c)
   h=freq/f_gyro
   wp=sqrt(4.*!dpi*Nel_v*e^2/m0) ; plasma freq in Hz
   omega= 2.*!dpi*freq           ; Freq. angular em ciclos/seg
   time0=systime(/sec)
   size=px/AU*180./!dpi*3600.    ; Size do pixel em arcsec

;   FOR i=0l,nn-1l DO BEGIN
;      Nel=Nel_v[i]
;      Np=Np_v[i]
;      B=Bmag_v[i]
;      ANGLE=ANGLE_v[i]
;      IF freq/(2.8e6*bmag_v) lt 500 THEN BEGIN
   plow=where(freq/(2.8e6*bmag_v) lt 30000,compl=phigh)   ; bom p/ b>10G e freq<1000GHz

;      GyroSyncs, Delta, Energy, Nel, M, B, Np, ANGLE, Freq, Flux, $
;                 Theta=Theta, LSource=LSource,SFU=1,phi1=phi1,	$
;		 phi2=phi2,Qstokes=Qstokes,Vstokes=Vstokes,	$
;		 e1d=e1d,e2d=e2d,a1d=a1d,a2d=a2d,SQUARE=1,	$
;		 TEMPERATURE=1E6,TRANSFRAD=1,/bessel,/quiet

if plow[0] ne -1 then begin
   bmag=bmag_v[plow]
   nel=nel_v[plow]
   angle=angle_v[plow]
   np=np_v[plow]
         gyroim_mp,freq,flux $
              ,delta=delta $
              ,energy=energy $
              ,nel=nel_v $ ; dens_NT $
              ,anor=anor $
              ,ntot=ntot $
              ,m=m $
              ,bmag=bmag $
              ,np=np $
              ,alpha=alpha $
              ,angle=angle $
              ,size=size $
              ,height=LSource $
              ,temperature=temperature $
              ,e1d=e1d,e2d=e2d,a1d=a1d,a2d=a2d $
              ,phi1=phi1, phi2=phi2 $
              ,electron=electron $
              ,plot=plot $ 
              ,oplot=oplot $
              ,QUIET=1 $
              ,ctime=ctime $ 
              ,single_core=single_core $  
              ,_EXTRA=_EXTRA 

         freq=freq[0] ; o que sera isso?
         eord[plow]=e1d
         extr[plow]=e2d
         aord[plow]=a1d
         axtr[plow]=a2d
ENDIF ELSE BEGIN
   if phigh[0] ne -1 then begin
         eord[phigh]=0.
         extr[phigh]=0.
         aord[phigh]=0.
         axtr[phigh]=0.
         flux=0.
   endif
endelse
;      ENDELSE 
      fluxo=fluxo*flux
;   ENDFOR ; esse e o for

   pp=where(fluxo eq 1.)
   if (pp ne -1) then fluxo[pp]=0.

   ress.B=Bmag_v
   ress.angle=angle_v
   ress.Nel=Nel_v
   ress.h=h
   ress.w_gyro=w_gyro
   ress.w_p=wp
   ress.eord=eord
   ress.extr=extr
   ress.aord=aord
   ress.axtr=axtr
;   ress.flux=fluxo
   res={gyro:in_gyro,res:ress}

end
