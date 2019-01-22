PRO read_gyro_im_mp, r, file=file

 IF ~keyword_set(file) THEN file='gyroim_mp.out'
 openr,unit,file,/get_lun
 IF ~eof(unit) THEN BEGIN 
    readf,unit,nfq
    readf,unit,fq
    jo=dblarr(nfq)
    jx=dblarr(nfq)
    ko=dblarr(nfq)
    kx=dblarr(nfq)

    i=0
    tmp=dblarr(4)
    WHILE ~eof(unit) DO BEGIN 
       readf,unit,tmp
;      fq[i]=tmp[0]
       jo[i]=tmp[0]
       jx[i]=tmp[1]
       ko[i]=tmp[2]
       kx[i]=tmp[3]
       i++
    ENDWHILE 
 ENDIF 
 free_lun,unit
; print,'freq=',fq,jo
 r = {fq:fq,jo:jo,ko:ko,jx:jx,kx:kx}
 
END 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO gyroim_mp,freq,flux $
          ,delta=delta $
          ,energy=energy $
          ,nel=nel $
          ,anor=anor $
          ,ntot=ntot $
          ,m=m $
          ,bmag=bmag $
          ,np=np $
          ,alpha=alpha $
          ,angle=angle $
          ,size=size $
          ,height=height $
          ,temperature=temperature $
          ,e1d=e1d,e2d=e2d,a1d=a1d,a2d=a2d $
          ,phi1=phi1, phi2=phi2 $
          ,electron=electron $
          ,plot=plot $ 
          ,oplot=oplot $
          ,QUIET=QUIET $
          ,ctime=ctime $ 
          ,single_core=single_core $
;          ,mode2=mode2 $  
          ,_EXTRA=_EXTRA

 if keyword_set(mode2) then fname=mode2 else fname=''

 ;; Constants:
 AU = 1.49597870d13                ; Astronomic Unit
 arc2cm = !dtor/3600d0 * AU        ; arcsec to cm in Sun
 m0 = 9.1094d-28                   ; electron mass [g]
 c = 2.998d10                      ; [cm/s] speed of light
 e = 4.803d-10                     ; electron charge
 E0 = m0 * c^2. / 1.6022d-12 / 1d3 ; electron rest energy [keV]

 if not keyword_set(delta) then delta=3.0
 if not keyword_set(energy) then energy=[10,5e3]
 if not keyword_set(m) then m=fix(0) ;; in radians
 if not keyword_set(size) then size=12.
 if not keyword_set(height) then height=5e8
 if not keyword_set(temperature) then temperature=0.0 ; t=0 e' o mesmo q.desligar o bremss

 nbmag=n_elements(bmag)

; vol = Height * !dpi*(size/2.*arc2cm)^2
 vol = Height * (size*arc2cm)^2

IF n_elements(delta) EQ 1 THEN delta=delta[0]

 IF ~keyword_set(Nel) AND ~keyword_set(Ntot) AND ~keyword_set(Anor) THEN $
  nel=1e7

IF n_elements(ntot) NE 0 THEN BEGIN  
    Nel  = Ntot / vol
    anor  = Ntot * (delta-1) / ((Energy[0]/1d3)^(-delta+1) - (Energy[1]/1d3)^(-delta+1)) 
 ENDIF ELSE BEGIN 
    IF n_elements(nel) NE 0 THEN BEGIN  
       Ntot  = Nel * vol
       anor  = Ntot * (delta-1) / ((Energy[0]/1d3)^(-delta+1) - (Energy[1]/1d3)^(-delta+1)) 
    ENDIF ELSE BEGIN 
       IF n_elements(anor) NE 0 THEN BEGIN 
          IF n_elements(energy) GT 2 THEN $
           message,'Anor input only for single power-law. Input Ntot or Nel instead.'
          ;; Anor is the number of electrons at 1 MeV
          ;; since we use energy in keV, we divide by 1E3
          Ntot = anor * ((Energy[0]/1d3)^(-delta+1) - (Energy[1]/1d3)^(-delta+1)) / (delta-1) 
          nel = Ntot / vol
       ENDIF 
    ENDELSE 
 ENDELSE 

; IF keyword_set(alpha) THEN BEGIN 
;    vb = 0.5 * e / (!pi * c * m0) * Bmag ; gyrofrequency
;    vp = 1.5*vb / alpha
;    np = vp^2 * !dpi *m0 / e^2
; ENDIF 

;; gyro.cpp stuff:

 ned = n_elements(energy)
 nfreq = n_elements(freq)
 if nfreq eq 0 or nfreq gt 1 then begin
    freq=1.7e10
    nfreq=1
 endif

 phi=0
 gphi=0


 IF n_elements(m) EQ 2 THEN BEGIN 
 ;;   print,'gaussian anisotropy'
    p0 = m[0]
    p1 = m[1]
    pitchangle=[p0,p1]
    mm = fix(1)
    npd = fix(2)
 ENDIF 
 IF n_elements(m) GT 2 THEN BEGIN 
  ;;  print,'array anisotropy'
    mm=fix(2)
    phi=m[*,0]
    gphi=m[*,1]
    npd=n_elements(phi)
    pitchangle=[[phi],[gphi]]
    ;; p0=!pi/2.*0. ;; temp
    ;; p1=!pi/6. ;; temp
    ;; npd = 100
    ;; phi=interpol([0.0,!pi],npd)
    ;; gphi=exp(-0.5*(phi-p0)^2/(p1)^2)
 ENDIF 
 IF n_elements(m) LE 1 THEN BEGIN
   ;; print,'isotropic'
    mm=0
    npd=2 ;; dummy values
 ENDIF 

; Garantia de angulo valido
 tmp_angle=angle    
 pangl=where(angle lt 3 and angle ge 0)
 if pangl[0] ne -1 then tmp_angle[pangl]=3.
 pangl=where(angle lt 0 and angle gt -3)
 if pangl[0] ne -1 then tmp_angle[pangl]=-3.
 pangl=where(angle gt 177. and angle le 180.)
 if pangl[0] ne -1 then tmp_angle[pangl]=177.
 pangl=where(angle gt 87. and angle le 90.)
 if pangl[0] ne -1 then tmp_angle[pangl]=87.

 if keyword_set(mode2) then openw,unit,'~/gyrotmp/'+fname+'gyroim.in',/get_lun else $
 openw,unit,'gyroim.in',/get_lun
 printf,unit, n_elements(bmag)
; printf,unit, nel
 printf,unit, [ned,mm,npd]
 IF nfreq NE 0 THEN printf,unit, [freq] 
 printf,unit, [energy, delta]
 printf,unit, bmag
 printf,unit, tmp_angle
 printf,unit, np
 printf,unit, nel
 IF mm NE 0 THEN printf,unit,pitchangle

 free_lun,unit

 time0=systime(/sec)

; diret="'/Users/usuario/Dropbox/ProgramasIDL/Tratamento_NORH/gyroim.out'"
; fim=" -i '/Users/usuario/Dropbox/ProgramasIDL/Tratamento_NORH/gyroim.in'"

 if keyword_set(mode2) then command = 'gyroim -o ~/gyrotmp/'+fname+'gyroim.out' else $
 command = 'gyroim_mp'; ///// -o '+diret+fim



 spawn, command, mensa ;nulls
;print,mensa
 ctime=systime(/sec)-time0

 IF ~keyword_set(QUIET) THEN print,'Elapsed time:',ctime,' seconds.'


 read_gyro_im_mp, r, file='gyroim_mp.out'

 freq=r.fq

 cgs2sfu = 1e19
 ;; circular area
 ; area = !pi * (size/2.0*arc2cm)^2.
 ;; square area
 area = (size*arc2cm)^2.
 volume = area * height
 ;; Source solid angle:
 omega = area / AU^2.
 ;; square area
 ;;th = size * !dtor / 3600.
 ;;omega = th^2.

 IF KEYWORD_SET(TEMPERATURE) THEN BEGIN
    ;; free-free coefficient by Dulk (1985).
    ;; assuming fully ionized hidrogen isothermal plasma.
    ;; ASSUMING HOMOGENEOUS DISTRIBUTION +++ NEED MODIFICATION!!!!!
    TEMPERATURE0=max(TEMPERATURE)
    kb = 1.38000e-016           ;; Boltzmann constant
    t1 = (TEMPERATURE0 LT 2e5) ? 18.2 + alog(TEMPERATURE0^(1.5)) - alog(freq) : $
         24.5 + alog(TEMPERATURE0) - alog(freq)
    kff = 9.78e-3 * np^2 / freq^2 / TEMPERATURE0^(1.5) * t1
    jff = kff * kb * TEMPERATURE0 * FREQ^2 / C^2
    r.jo += jff
    r.jx += jff
    r.ko += kff
    r.kx += kff
;;  print,'FREE-FREE EMISSION/ABSORPTION INCLUDED.'
 ENDIF

;; RADIATIVE TRANSFER
 phi1=dblarr(n_elements(freq))
 phi2=dblarr(n_elements(freq))
 tau_o=r.ko*height
 tau_x=r.kx*height

 toLow = where(tau_o LT 0.001)
 toMed = where((tau_o GE 0.001) AND (tau_o LE 20.))
 toHig = where(tau_o GT 20.)
 txLow = where(tau_x LT 0.001)
 txMed = where((tau_x GE 0.001) AND (tau_x LE 20.))
 txHig = where(tau_x GT 20.)
 IF toLow[0] NE -1 THEN phi1[toLow] = r.jo[toLow] * omega * height
 IF toMed[0] NE -1 THEN phi1[toMed] = omega * r.jo[toMed]/r.ko[toMed]*(1d0-exp(-tau_o[toMed]))
 IF toHig[0] NE -1 THEN phi1[toHig] = r.jo[toHig]/r.ko[toHig] * omega
 IF txLow[0] NE -1 THEN phi2[txLow] = r.jx[txLow] * omega * height
 IF txMed[0] NE -1 THEN phi2[txMed] = omega * r.jx[txMed]/r.kx[txMed]*(1d0-exp(-tau_x[txMed]))
 IF txHig[0] NE -1 THEN phi2[txHig] = r.jx[txHig]/r.kx[txHig] * omega
 
;; phi1 = r.jo/r.ko*(1.0-exp(-r.ko*height))*omega*cgs2sfu
;; phi2 = r.jx/r.kx*(1.0-exp(-r.kx*height))*omega*cgs2sfu
 phi1*=cgs2sfu
 phi2*=cgs2sfu
 flux = phi1+phi2

 ;; phi1 = r.jo/r.ko*(1.0-exp(-r.ko*height))*omega*cgs2sfu
 ;; phi2 = r.jx/r.kx*(1.0-exp(-r.kx*height))*omega*cgs2sfu
 ;; flux = phi1+phi2

 e1d=r.jo
 e2d=r.jx
 a1d=r.ko
 a2d=r.kx

 electron={nel:nel,ntot:ntot,delta:delta,energy:energy}
 
 if keyword_set(plot) then BEGIN
    yr=[min(flux[where(flux GT 0)]),max(flux[where(flux LT 1e10)])]
    if n_elements(freq) lt 10 then psym=-4
    plot,freq/1e9,flux,/xl,/yl,psym=psym ,_EXTRA=_EXTRA $
         ,xtit='Frequency [GHz]',ytit='Flux Density [sfu]',yr=yr
    if n_elements(flux_old) gt 0 THEN $
     oplot,freq/1e9,flux_old,col=fsc_color('red'),lines=2 ,_EXTRA=_EXTRA
 ENDIF

 IF keyword_set(oplot) THEN oplot,freq,flux,_EXTRA=_EXTRA

END 
