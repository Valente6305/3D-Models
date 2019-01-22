;+
; NAME:search_catalogue
;
;
;
; PURPOSE: Find the best fit spectra in the catalogue for the 
;          flux density given (as input)
;
;
;
; CATEGORY:
;
;
;
; CALLING SEQUENCE:search_catalogue,sfu,hel [,deltalat=deltalat,deltalong=deltalong]
;
;
;
; INPUTS: sfu = a vector with 4 elements (the observed flux density at
;               3.75,9.4,17,34 GHz)
;         hel = a string with the burst (active region) coordinate,
;         ex: 'N07W32'
;
;
;
; OPTIONAL INPUTS: if keyword_set(print1)  azimuth equal (az=-60),print2 (az=-15),
; print3 (15),print4 (az=60)
;
;
;
; KEYWORD PARAMETERS:deltalat and deltalong = define the sides of a box in latitude
; and longitude (deg) to serach for models. Default numbers are 10 deg
; x 10 deg
;
;
;
; OUTPUTS: page in the calalogue of the best model to fit data
;
;
;
; OPTIONAL OUTPUTS:
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
; EXAMPLE: search_catalogue,[71.,255.,386.,172],'S18E33'
;         [Example: Nobeyama Radioheliograph Event 20120710_0339]
;
;
;
; MODIFICATION HISTORY:
; J.E.R Costa, Nov,2012
;V.A. Cuambe, Nov,2014
;V.A. Cuambe, Abril,2015
pro search_catalogue,sfu,hel,res,deltalat=deltalat,deltalong=deltalong,print0=print0,print1=print1,print2=print2,print3=print3,print4=print4,plot=plot,catal=catal,n_mod=n_mod,ii=ii
if total(size(catal)) lt 100 then restore,'Catalogue_completo.sav' else st_models=catal
if total(size(n_mod)) eq 0 then restore,'Catalogue_completo.sav'
catal=st_models
n=n_elements(st_models)
if not keyword_set(deltalat) then dlat=10.0 else dlat=deltalat
if not keyword_set(deltalong) then dlon=10.0 else dlon=deltalong
if strpos(hel,'S') ne -1 then slat=-1. else slat=1.0
if strpos(hel,'E') ne -1 then slon=-1. else slon=1.0

b_lat=slat*float(strmid(hel,1,2))
phi=slon*float(strmid(hel,4,2))

coord=[phi,b_lat]
freq=[3.75,9.4,17,34 ]

pp=where(abs(st_models.phi-coord[0]) lt dlon and abs(st_models.b_lat-coord[1]) le dlat)
nn=n_elements(pp)
dd=total( (alog10(st_models[pp].flux)-alog10(sfu#replicate(1.,nn)))^2,1)
ss=sort(dd)
lim=(n_elements(dd)-1)<11
sigma=dd[ss[0:lim]]
flux=st_models[pp[ss[0:lim]]].flux
res={freq:freq,flux:flux,sigma:sigma}
if n_elements(sfu) ne 4 then begin
   print,'Please enter the flux densities of the [3.75,9.4,17,34 ] frequencies!'
   return
endif
cor=[30,50,70,95,110,130,150,180,200,220,240,250]
if keyword_set(plot) then begin
; Finding models near location
   window,0,xsiz=300,ysiz=400
   !p.multi=0   ;[0,4,3]
   for k=0,11 do begin 
      if k eq 0 then plot,freq,sfu,psy=2,/xlog,/ylog,xtit='Frequency (GHz)',ytit='Flux density(sfu)'
      oplot,freq,st_models[pp[ss[k]]].flux,col=cor[k]
      ye=10.^!y.crange
      py=10.^((!y.crange[1]-!y.crange[0])*0.95+!y.crange[0]-(!y.crange[1]-!y.crange[0])/35.*k)
      oplot,[1.4,2.0],[py,py],col=cor[k]
      xyouts,1.2,py,string(k+1,format='(i2)'),alignm=0.5
   endfor     
   window,1,tit='17GHz',xsiz=4l*64,ysiz=3l*64
   for k=0,11 do begin 
      www=bytscl(rebin(st_models[pp[ss[k]]].image[*,*,2],64,64))
      www[0,*]=255
      www[63,*]=255
      www[*,0]=255
      www[*,63]=255
      tv,www,k ; 17 GHz
   endfor     
endif
;plot,freq,sfu,psy=2,/xlog,/ylog
;oplot,freq,st_models[pp[ss[1]]].flux      
;plot,freq,sfu,psy=2,/xlog,/ylog
;oplot,freq,st_models[pp[ss[2]]].flux 
tt=min((pp[ss]-total(n_mod,/cum)) gt 0,ichap)
q=where(ichap eq 0)
if q[0] ne -1 then ichap[q]=2
ichap=ichap-1     
pages=pp[ss]-total(n_mod,/cum)
ii=pp[ss[0:lim]]
!p.multi=0
fmt='(i2,a1,2x,i6,a1,2x,a4,f5.1,2x,a4,f5.1,2x,a4,f5.1,2x,a5,f5.1,2x,a5,f5.1,2x,f6.4)'
if keyword_set(print0) then begin
   print,'Twelve possibilities are given:>'
   print,' #     index    long        lat      azimuth    Inclin    Assymet     Chi2'

   for i=0,20 do begin 
      print,i+1,': ',pp[ss[i]],':','Phi:',St_models[pp[ss[i]]].phi,'Lat:',st_models[pp[ss[i]]].b_lat,'Az:',st_models[pp[ss[i]]].az,'Incl:',st_models[pp[ss[i]]].in, 'Ass: ',st_models[pp[ss[i]]].as, dd[ss[i]],format=fmt
   endfor
endif
if keyword_set(print2) then begin
   print,'Twelve possibilities are given at preferred azimuth:>'
   print,' #     index    long        lat      azimuth    Inclin    Assymet     Chi2'
      
   for i=0,20 do begin
 j=where(st_models[pp[ss]].az eq -15)
     print,j[i]+1,': ',pp[ss[j[i]]],':','Phi:',St_models[pp[ss[j[i]]]].phi,'Lat:',st_models[pp[ss[j[i]]]].b_lat,'Az:',st_models[pp[ss[j[i]]]].az,'Incl:',st_models[pp[ss[j[i]]]].in, 'Ass: ',st_models[pp[ss[i]]].as,dd[ss[j[i]]],format=fmt
  endfor
endif
if keyword_set(print3) then begin
   print,'Twelve possibilities are given at preferred azimuth:>'
   print,' #     index    long        lat      azimuth    Inclin    Assymet     Chi2'
   
   for i=0,20 do begin
 j=where(st_models[pp[ss]].az eq 15)
     print,j[i]+1,': ',pp[ss[j[i]]],':','Phi:',St_models[pp[ss[j[i]]]].phi,'Lat:',st_models[pp[ss[j[i]]]].b_lat,'Az:',st_models[pp[ss[j[i]]]].az,'Incl:',st_models[pp[ss[j[i]]]].in, 'Ass: ',st_models[pp[ss[i]]].as,dd[ss[j[i]]],format=fmt
  endfor
endif
if keyword_set(print1) then begin
   print,'Twelve possibilities are given at preferred azimuth:>'
   print,' #     index    long        lat      azimuth    Inclin    Assymet     Chi2'
   print,''
   for i=0,20 do begin
 j=where(st_models[pp[ss]].az eq -60)
     print,j[i]+1,': ',pp[ss[j[i]]],':','Phi:',St_models[pp[ss[j[i]]]].phi,'Lat:',st_models[pp[ss[j[i]]]].b_lat,'Az:',st_models[pp[ss[j[i]]]].az,'Incl:',st_models[pp[ss[j[i]]]].in, 'Ass: ',st_models[pp[ss[i]]].as,dd[ss[j[i]]],format=fmt
  endfor
endif
if keyword_set(print4) then begin
   print,'Twelve possibilities are given at preferred azimuth:>'
   print,' #     index    long        lat      azimuth    Inclin    Assymet     Chi2'

   for i=0,20 do begin
 j=where(st_models[pp[ss]].az eq 60)
     print,j[i]+1,': ',pp[ss[j[i]]],':','Phi:',St_models[pp[ss[j[i]]]].phi,'Lat:',st_models[pp[ss[j[i]]]].b_lat,'Az:',st_models[pp[ss[j[i]]]].az,'Incl:',st_models[pp[ss[j[i]]]].in, 'Ass: ',st_models[pp[ss[i]]].as,dd[ss[j[i]]],format=fmt
  endfor

endif
end
