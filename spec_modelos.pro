pro spec_modelos
 Restore,'flux_distribution4.sav'
 close,1
 ii=0l
 in=0l
 openw,1,'all_thick.txt'
 n=n_elements(st)
 freq=[2.5,3.70,9.30,17.00,35.0]
 w=fltarr(5)+1.0
 w=[1.,30.,50.,20.0,1.0]
 a=fltarr(n,5)
 fmax=fltarr(n)
 
 for i=0l, n-1l do begin
    f=[0.,reform(flux[*,i])]
    alfa=alog(f[2]/f[1])/alog(9.3/3.7)

    f[0]=f[1]*((2.5/3.7)^alfa)
    yy=linfit(alog(freq[1:2]),alog(f[1:2]))
    zz=linfit(alog(freq[3:4]),alog(f[3:4]))
    vp=(zz[0]-yy[0])/(yy[1]-zz[1])
    at=[1.1*exp(zz[1]*vp+zz[0]),1.1*zz[1],exp(vp),1.1*alfa]
;print,at
    result=spec_fit(freq,f,w,at,$
          function_name='espec5',iter=niter,itmax=100,tol=0.01,flag=flag,chi2=chi2)
    a[i,*]=[at,chi2]
    x=findgen(100)/99.*48.+2.
    espec5,x,at,y
    tt=max(y,iii)
    fmax[i]=x[iii]

;print,[at,sqrt(chi2/5.)/at[0]]
;    Res = POLY_FIT( freq,reform(flux[*,i]) , 2)

goto,pula
if sqrt(chi2/5.)/at[0] lt 1000.8 then begin
    plot,freq,f,/xlog,/ylog,psy=2
    x=findgen(30)/29.*48.+2.
;    y=res[0]+res[1]*x+res[2]*x^2 ;+res[3]*x^3
    espec5,x,at,y
    oplot,x,y
 ;   stop
wait,0.5
endif
pula:
if at[1] gt 0.0 then begin
   in=[in,i]

;   plot,freq,f,/xlog,/ylog,psy=2
;   x=findgen(30)/29.*48.+2.
;   espec5,x,at,y
;   oplot,x,y
;   wait,0.2

   printf,1,ff[i]
endif else begin 

;   plot,freq,f,/xlog,/ylog,psy=2
;   x=findgen(30)/29.*48.+2.
;   espec5,x,at,y
;   oplot,x,y
;   wait,0.2

   ii=[ii,i]
endelse
 endfor
close,1
ii=ii[1:*]
in=in[1:*]

stthick=st[in]
ffthick=ff[in]
stpeak=st[ii]
ffpeak=ff[ii]
save,file='resthick4.sav',ii,in,stpeak,ffpeak,stthick,ffthick,a,fmax
end
