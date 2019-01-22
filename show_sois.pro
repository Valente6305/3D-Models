pro show_sois,nome,ps=ps  ;,fig,r,g,b
restore,nome
DC=0.22
mydevice=!D.NAME
loadct,3,/silent
tvlct,r,g,b,/get
st=mag_res.st_arc
show_sol,st,fig,r,g,b
rad='/Users/usuario/Dropbox/ProgramasIDL/Tratamento_NORH/Figuras/'
p=strpos(nome,'nn')
q=strpos(nome,'.sav')
nn=rad+'sol'+strmid(nome,p+2,q-p-2)+'.png'
if keyword_set(ps) then begin
   mydevice = !D.NAME
   ; Set plotting to PostScript:
   SET_PLOT, 'PS'
   tvlct,r,g,b
   tv,fig[*,85:315],0.0231,0.627+DC,/norm,xsize=0.16,ysize=.17
endif else begin
   write_png,nn,fig,r,g,b
endelse
show_sol,st,fig,r,g,b,/z
nn=rad+'loop'+strmid(nome,p+2,q-p-2)+'.png'
if keyword_set(ps) then begin
   mydevice = !D.NAME
   ; Set plotting to PostScript:
   SET_PLOT, 'PS'
   tvlct,r,g,b
   tv,fig[*,20:380],0.0351,0.40+DC,/norm,xsize=0.14,ysize=.195
endif else begin
   write_png,nn,fig,r,g,b
endelse
set_plot,mydevice
end
