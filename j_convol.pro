pro j_convol,imagem,w17,w35
   AU = 1.49597870d13 
   arc2cm = !dtor/3600d0 * AU        ; arcsec to cm in Sun
   pix=imagem.gyro.px
   arc=15.0
   dres17=arc*arc2cm/pix/sqrt(alog(256.))   ; resolucao angular em pixels
   arc=7.5
   dres35=arc*arc2cm/pix/sqrt(alog(256.))
   w17=total(imagem.image[2].imagem,1)
   w35=total(imagem.image[3].imagem,1)
   wid=29
   kernel17=exp(-(shift(dist(wid),wid/2,wid/2))^2/dres17^2/2.)
   kernel35=exp(-(shift(dist(wid),wid/2,wid/2))^2/dres35^2/2.)
   w17=convol(w17,kernel17,/edge_truncate)
   w35=convol(w35,kernel35,/edge_truncate)
;   w17=gauss_smooth(w17,dres17,width=29,/edge_truncate,kernel=v)
;   w35=gauss_smooth(w35,dres35,width=29,/edge_truncate,kernel=v)

;   CONTOUR, rebin(reform(w[*,*,ii]),256,256), PATH_XY=xy, PATH_INFO=info ,levels=[max(w[*,*,ii])/2.]

end
