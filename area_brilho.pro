pro area_brilho, img,flux,ii,px=px,res=res,area17=area17,conv=conv,f17=f17
   if not keyword_set(res) then arc=15.0 else arc=res
   if not keyword_set(px) then pix=7.34116e+07 else pix=px
   
   AU = 1.49597870d13 
   arc2cm = !dtor/3600d0 * AU        ; arcsec to cm in Sun

   dres=arc*arc2cm/pix   ; resolucao angular em pixels
;   print,'HPBW(PIXELS)=',dres
   dres=dres/sqrt(alog(256.))  ; resolucao em sigma e nao HPBW (mult.p/2.35)
;   print,'Sigma=',dres
   if total(size(ii)) eq 0 then ii=2   ; era 3??!
   flux=fltarr(4)
   xx=findgen(256)*pix/8./arc2cm
   w=total(img.image.imagem,1)
   if keyword_set(conv) then begin 
      wid=16
      kernel=exp(-(shift(dist(wid),wid/2,wid/2))^2/dres^2/2.)
      w[*,*,ii]=convol(reform(w[*,*,ii]),kernel,/edge_truncate)
;      w[*,*,ii]=gauss_smooth(w[*,*,ii],dres,width=29,/edge_truncate,kernel=v)
      if ii eq 2 then f17=total(w[*,*,2])
   endif
;stop
   denovo:
   CONTOUR, rebin(reform(w[*,*,ii]),256,256),xx,xx, PATH_XY=xy, /PATH_DATA_COORDS, PATH_INFO=info ,levels=[max(w[*,*,ii])/2.]
   n=n_elements(info)
   if n eq 0 then begin
      print,'Problema de NAN na imagem da freq[] ',ii
      qqq=reform(w[*,*,ii])
      pppp=where(finite(qqq) eq 0)
      qqq[pppp]=0.0
      w[*,*,ii]=qqq
      goto, denovo
   endif
   area=fltarr(n)
   in=total(info.n,/cum)-info.n
   fim=total(info.n,/cum)-1l
;   if total(info.high_low) ne n then stop
   goto,pula
 ;  w[*,*,0]=gauss_smooth(w[*,*,0],dres,width=29,/edge_truncate,kernel=v)
 ;  w[*,*,1]=gauss_smooth(w[*,*,1],dres,width=29,/edge_truncate,kernel=v)
  ; w[*,*,2]=gauss_smooth(w[*,*,2],dres,width=29,/edge_truncate,kernel=v)
  ; w[*,*,3]=gauss_smooth(w[*,*,3],dres,width=29,/edge_truncate,kernel=v)
;   goto, pula
   !p.multi=[0,1,4]
   contour,rebin(reform(w[*,*,0]),256,256),/iso,levels=[max(w[*,*,0])/2.]
   contour,rebin(reform(w[*,*,1]),256,256),/iso,levels=[max(w[*,*,1])/2.]
   contour,rebin(reform(w[*,*,2]),256,256),/iso,levels=[max(w[*,*,2])/2.]
   contour,rebin(reform(w[*,*,3]),256,256),/iso,levels=[max(w[*,*,3])/2.]

   pula:
  
   for i=0,n-1l do begin
      area[i]=POLY_AREA(XY[0,in[i]:fim[i]], XY[1,in[i]:fim[i]])
   endfor

;   print,'areas=',area
;   print,'total=',total(area)
   flux=total(total(w,1),1)
   
   area17=total(area)

end
