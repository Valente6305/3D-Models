PRO read_gyro, r, file=file

 IF ~keyword_set(file) THEN file='gyro.out'
 openr,unit,file,/get_lun
 IF ~eof(unit) THEN BEGIN 
    readf,unit,nfq
    
    fq=dblarr(nfq)
    jo=dblarr(nfq)
    jx=dblarr(nfq)
    ko=dblarr(nfq)
    kx=dblarr(nfq)

    i=0
    tmp=dblarr(5)
    WHILE ~eof(unit) DO BEGIN 
       readf,unit,tmp
       fq[i]=tmp[0]
       jo[i]=tmp[1]
       jx[i]=tmp[2]
       ko[i]=tmp[3]
       kx[i]=tmp[4]
       i++
    ENDWHILE 
 ENDIF 
 free_lun,unit
 
 r = {fq:fq,jo:jo,ko:ko,jx:jx,kx:kx}
 
END 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO gyro,freq,flux $
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
          ,keep=keep $
          ,new=new $
          ,_EXTRA=_EXTRA

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
 if not keyword_set(bmag) then bmag=500.
 if not keyword_set(angle) then angle=45.
 if not keyword_set(np) then np=1e9
 if not keyword_set(size) then size=12.
 if not keyword_set(height) then height=5d8
 if not keyword_set(temperature) then temperature=0.0


 vol = Height * !dpi*(size/2.*arc2cm)^2

IF n_elements(delta) EQ 1 THEN delta=delta[0]

 IF ~keyword_set(Nel) AND ~keyword_set(Ntot) AND ~keyword_set(Anor) THEN $
  nel=1d7

IF n_elements(ntot) NE 0 THEN BEGIN 
   ntot=double(ntot)
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

 IF keyword_set(alpha) THEN BEGIN 
    vb = 0.5 * e / (!pi * c * m0) * Bmag ; gyrofrequency
    vp = 1.5*vb / alpha
    np = vp^2 * !dpi *m0 / e^2
 ENDIF 

;; gyro.cpp stuff:

 ned = n_elements(energy)
 nfreq = n_elements(freq)

 phi=0
 gphi=0

 IF n_elements(m) EQ 2 THEN BEGIN 
    ;;print,'gaussian anisotropy'
    p0 = m[0]
    p1 = m[1]
    pitchangle=[p0,p1]
    mm = fix(1)
    npd = fix(2)
 ENDIF 
 IF n_elements(m) GT 2 THEN BEGIN 
    ;;print,'array anisotropy'
    size_m = size(m)
    IF size_m[2] EQ 2 THEN BEGIN 
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
    ENDIF ELSE BEGIN 
       mm=fix(4)
       npd=fix(m[0])
       pitchangle=m[1:*]
    ENDELSE 
 ENDIF 
 IF n_elements(m) LE 1 THEN BEGIN
    ;;print,'isotropic'
    mm=0
    npd=2 ;; dummy values
 ENDIF 
   
 openw,unit,'gyro.in',/get_lun
 printf,unit, [bmag,angle,np,nel]
 printf,unit, [ned,mm,npd,nfreq]
 IF nfreq NE 0 THEN printf,unit, freq 
 printf,unit, [energy, delta]
 IF mm NE 0 THEN printf,unit,pitchangle
 
 free_lun,unit

 time0=systime(/sec)

 ;command = '/home/paulo/solar/gyro/gyro'
 cd,cur=adir
 command='gyro -i '+adir+'/gyro.in -o '+adir+'/gyro.out'
 spawn, command, nulls, exit_status=exit_status,stderr=stderr
 ctime=systime(/sec)-time0

 IF ~keyword_set(QUIET) THEN print,'Elapsed time:',ctime,' seconds.'

 read_gyro, r

 ;; remove input and output ascii files after running the program (set
 ;; /keep to keep the files).
 if ~keyword_set(keep) then file_delete,'gyro.in','gyro.out',/ALLOW_NONEXISTENT,/NOEXPAND_PATH,/quiet

 freq=r.fq

 cgs2sfu = 1e19
 ;; circular area
 area = !pi * (size/2.0*arc2cm)^2.
 volume = area * height
 ;; Source solid angle:
 omega = area / AU^2.
 ;; square area
 ;;th = size * !dtor / 3600.
 ;;omega = th^2.

 IF KEYWORD_SET(TEMPERATURE) THEN BEGIN
    ;; free-free coefficient by Dulk (1985).
    ;; assuming fully ionized hidrogen isothermal plasma.
    kb = 1.38000e-016           ;; Boltzmann constant
    t1 = (TEMPERATURE LT 2e5) ? 18.2 + alog(TEMPERATURE^(1.5)) - alog(freq) : $
         24.5 + alog(TEMPERATURE) - alog(freq)
    kff = 9.78e-3 * np^2 / freq^2 / TEMPERATURE^(1.5) * t1
    jff = kff * kb * TEMPERATURE * FREQ^2 / C^2
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
    if not keyword_set(new) then $
       plot,freq/1e9,flux,/xl,/yl,psym=psym ,_EXTRA=_EXTRA $
         ,xtit='Frequency [GHz]',ytit='Flux Density [sfu]',yr=yr  ; $
;    else p=plot(freq/1e9,flux,/xlog,/ylog,'-r2' ,_EXTRA=_EXTRA $
;         ,xtitle='Frequency [GHz]',ytitle='Flux Density [sfu]',font_style='bf',yrange=[1.0,yr[1]];)
    if n_elements(flux_old) gt 0 THEN $
     oplot,freq/1e9,flux_old,col=fsc_color('red'),lines=2 ,_EXTRA=_EXTRA
 ENDIF

 IF keyword_set(oplot) THEN begin 
    if not keyword_set(new) then oplot,freq/1e9,flux,_EXTRA=_EXTRA ;else $
;       q=plot(freq/1e9,flux,_EXTRA=_EXTRA,'-b2',/overplot)
 endif

END 
