function convol_17,imagem,px,resol=resol
   AU = 1.49597870d13 
   arc2cm = !dtor/3600d0 * AU        ; arcsec to cm in Sun
   pix=px ;imagem.gyro.px
   if not keyword_set(resol) then arc=12.0 else arc=resol
   dres17=arc*arc2cm/pix/sqrt(alog(256.))   ; resolucao angular em pixels
   w17=imagem
   wid=n_elements(imagem[*,0])/4*2+1   ;29
   kernel17=exp(-(shift(dist(wid),wid/2,wid/2))^2/dres17^2/2.)
   return,convol(w17,kernel17,/edge_truncate)
end
