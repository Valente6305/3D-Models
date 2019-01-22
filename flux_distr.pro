pro flux_distr,st,names,lat,long,flux,area,im,file=file
; Input st and names [st and names are outputs of file_distrib]
; output lat,long,flux, area
pp=where(st.freq eq '0037009301700350')
tnames=names[pp]
tst=st[pp]
n=n_elements(pp)
flux=fltarr(4,n)
ang=fltarr(n)
area=fltarr(6,n)  ; aumentei para conter a area de 17GHz conv c/ 15 arcsec
im=fltarr(32,32,4,n)
lat=tst.b_lat
long=tst.phi
tempo0=systime(/sec)
for i=0l, n-1l do begin
   restore,tnames[i] 
   for j=0,3 do begin  ; n# de freqs
      area_brilho,imagem,tflux,j,px=imagem.gyro.px,area17=a ; o j=2 é o 17 GHz
      area[j,i]=a
   endfor
   area_brilho,imagem,tflux,2,px=imagem.gyro.px,area17=a,/conv,res=15.0 ; o j=2 é o 17 GHz
   area[4,i]=a
   area_brilho,imagem,tflux,3,px=imagem.gyro.px,area17=a,/conv,res=8.0 ; o j=2 é o 17 GHz
   area[5,i]=a
   w=total(imagem.image.imagem,1)
   im[*,*,*,i]=w
   flux[*,i]=total(total(w,1),1)
   if i mod 1000 eq 0 then print,i,systime()
endfor
ff=names
save,file='flux_distribution.sav',st,ff,lat,long,flux,area,w  
if keyword_set(file) then save,file=file,st,ff,lat,long,flux,area,w  
end
