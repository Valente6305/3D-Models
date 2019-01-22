pro roda_imag,nome,ps=ps
restore,nome
if n_elements(nome) eq 0 then nome='Spectrum'
n=n_elements(imagem.image)
fluxos=fltarr(n)
nsq=ceil(sqrt(n))
spix=128
;window, xs=300, ys=400,/free,title=string(imagem[i].freq/1.e9)+' GHz'
fmt='('+strcompress(string(n)+'f5.1)',/remove_all)
tit=string(imagem.freq/1.e9,format=fmt)+' GHz'
;window, xs=spix*nsq, ys=spix*nsq,/free,title=tit 
set_plot,'Z'
DEVICE,set_RES=[300,400]      
tvlct,RRR,GGG,BBB,/get
peso=reform(imagem.image[*].imagem[0,*,*])*0.0
for i=0, n-1 do begin 
   ii=i/nsq
   jj=i mod nsq
   stokes_I=imagem.image[i].imagem[0,*,*]+imagem.image[i].imagem[1,*,*]
   stokes_V=imagem.image[i].imagem[0,*,*]-imagem.image[i].imagem[1,*,*]
   peso[*,*,i]=stokes_i
   im=bytscl(rebin(reform(stokes_I),spix,spix))
;   tv,im,jj*spix,(nsq-ii-1)*spix 
   fluxos[i]=total(stokes_I)
;   xyouts,jj*spix+fix(0.39*spix),(nsq-ii-1)*spix+fix(0.937*spix),fluxos[i],/dev
endfor

ttt=max(fluxos,i_max)
;=========
pesot=reform(peso[*,*,i_max]/total(peso[*,*,i_max]))
pesov=fltarr(32,32,32)
for j=0,31 do pesov[*,*,j]=pesot
depth=total(pesot*TOTAL((mag_res.b NE 0),3))*imagem.gyro.px
pp=where(mag_res.b ne 0)
pt=where(mag_res.theta gt 90)
if pt[0] ne -1 then mag_res.theta[pt]=180.-mag_res.theta[pt]
bm=total(mag_res.b[pp]*pesov[pp])/total(pesov[pp])
nel=total(mag_res.nel[pp]*pesov[pp])/total(pesov[pp])
angle=total(mag_res.theta[pp]*pesov[pp])/total(pesov[pp])

;=======
;pesot=reform(peso[*,*,i_max]/total(peso[*,*,i_max]))
;depth=total(pesot*TOTAL((mag_res.b NE 0),3))*imagem.gyro.px

;window, xs=spix*nsq, ys=spix*nsq,/free,title=nome
if keyword_set(ps) then set_plot,'PS'
;ff=findgen(100)/99.*88.+2
ff=10.^(findgen(50)/49.*1.6+0.3)
;pp=where(imagem.res.b ne 0)  ; Esta errado e tem que usar o mag_re.b e o peso
;bm=mean(imagem[pp].res.b) ; Esta errado e tem que usar o mag_re.b e o peso
;angle=mean(imagem[pp].res.angle)  ; Esta errado e tem que usar o mag_re.b e o peso
;nel=mean(imagem[pp].res.nel); Esta errado e tem que usar o mag_re.b e o peso

gyro,ff*1.e9,flux,delta=imagem.gyro.delta[0],angle=angle, size=10.0, height=depth, Nel=nel, $
     Bmag=bm,energy=imagem.gyro.energy,/quiet
tm=max(flux)
ttm=max(fluxos)
plot,imagem.freq/1.e9,fluxos,/xlog,/ylog,ytit='Flux Density (sfu)',xtit='Frequency (GHz)',psym=2, $
backg=255,col=0,position=[0.26,0.72,0.43,0.995],/noer,chars=0.5,title='D_factor='+ $
string(sqrt(ttm/tm),format='(f7.3)'),syms=0.5
oplot,ff,flux*ttm/tm,col=0 ;,pos=[0.65,.05,.95,.5]
q=[imagem.gyro.delta[0],angle, 10.0, depth, Nel,bm,imagem.gyro.energy]
;save,file='partmp.sav',q
;print,flux
if angle gt 90 then angle=180.-angle
label = string(bm,format='("B=",f6.1," G")')         + '!C!C' $
	      + string(ANGLE,format='("!4X!3 = ",f5.1," !eo!n")') + '!C!C' $
              + string(imagem.gyro.delta[0],format='("!4d!3=",f4.1)')          + '!C!C' $
	      + string(depth,format='("L=",e7.1," cm")')             + '!C!C' $
	      + string(Nel, format='("Nel=",e7.1, " cm!e-3!n")')

xyouts,.34,.83,label,col=0,/norm,charsiz=0.6

rad='~/Dropbox/idl_2015/';'/Users/usuario/Dropbox/ProgramasIDL/Tratamento_NORH/Figuras/'
p=strpos(nome,'nn')
q=strpos(nome,'.sav')
nn=rad+'spec'+strmid(nome,p+2,q-p-2)+'.png'
if not keyword_set(ps) then begin
   fig=tvrd()
   write_png,nn,fig,r,g,b
   set_plot,'X'
endif


end
