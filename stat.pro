pro stat
; Qual e o view primordial em cada freq?
; Qual a distr. de params nessa amostra?
; Pitch angle ditr. influencia no espec?

ff=findfile('Modelos2/img*.sav')
nf=n_elements(ff)
area=fltarr(6)
for i=0,nf-1 do begin
   restore,ff[i]
   w=double(reform(imagem.image.imagem[0,*,*]))
   fluxos=total(total(w,1),1)
   for j=0,n_elements(w[0,0,*])-1 do begin
      lev=0.5d0*max(w[*,*,j])
      contour,w[*,*,j],levels=[lev],path_xy=xy,/path_data_coord,/path_double,path_info=var
      tt=max(var.n,ii)
      if ii ne 0 then tot=total(var[0:ii-1].n) else tot=0
      xy=xy[*,tot:var[ii].n-1]
;      print,'j=',j, '-->Level=',lev,'info_lev=',var.level
      area[j]=poly_area(xy[0,*],xy[1,*])
   endfor
   plot,imagem.freq/1.e9,area,/xlog,/ylog, $
             ytit='Area (pixel 2)',xtit='Freq (GHz)',psym=2
wait,0.5
;stop
endfor

end
