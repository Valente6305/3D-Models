pro show_imagsp, imagem, nome=nome
if keyword_set(nome) then begin
   restore,nome
   pp=strpos(nome,'phi')
   nome=strmid(nome,pp,255)
endif
if n_elements(nome) eq 0 then nome='Spectrum'
n=n_elements(imagem.image)
fluxos=fltarr(n)
nsq=ceil(sqrt(n))
spix=128
;window, xs=300, ys=400,/free,title=string(imagem[i].freq/1.e9)+' GHz'
fmt='('+strcompress(string(n)+'f5.1)',/remove_all)
tit=string(imagem.freq/1.e9,format=fmt)+' GHz'
window, xs=spix*nsq, ys=spix*nsq,/free,title=tit 
for i=0, n-1 do begin 
   ii=i/nsq
   jj=i mod nsq
   stokes_I=imagem.image[i].imagem[0,*,*]+imagem.image[i].imagem[1,*,*]
   stokes_V=imagem.image[i].imagem[0,*,*]-imagem.image[i].imagem[1,*,*]
   im=bytscl(rebin(reform(stokes_I),spix,spix))
   tv,im,jj*spix,(nsq-ii-1)*spix 
   fluxos[i]=total(stokes_I)
   xyouts,jj*spix+fix(0.39*spix),(nsq-ii-1)*spix+fix(0.937*spix),fluxos[i],/dev
endfor

window, xs=spix*nsq, ys=spix*nsq,/free,title=nome
plot,imagem.freq/1.e9,fluxos,/xlog,/ylog,ytit='Flux Dens (sfu)',xtit='Freq (GHz)',psym=2
ff=findgen(100)/99.*88.+2
pp=where(imagem.res.b ne 0)
bm=mean(imagem[pp].res.b)
angle=mean(imagem[pp].res.angle)
nel=mean(imagem[pp].res.nel)
gyro,ff*1.e9,flux,delta=imagem.gyro.delta[0],angle=angle, size=10.0, height=imagem.gyro.px, Nel=nel, Bmag=bmag,/quiet
tm=max(flux)
ttm=max(fluxos)
oplot,ff,flux*ttm/tm ;,pos=[0.65,.05,.95,.5]

end
