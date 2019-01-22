pro show_loops,nome,fig,r,g,b,ps=ps
DC=0.22
mydevice=!D.NAME
print, 'Device=',mydevice
restore,nome
;window,xsize=800,ysize=400
nf=(size(imagem.image))[1]
im=bytarr(8,250,250)+255
im0=bytarr(250,250)
j_convol,imagem,w17,w35
;loadct,39  ;  AQUI

for i=0, nf-1 do begin 
   set_plot,'z'                  
   DEVICE,set_RES=[250,250]      
   if i eq 0 or i eq 1 then show_loop,nome=nome,ind=i else begin
      ww=w17
      if i eq 3 then ww=w35
      show_loop,nome=nome,ind=i,contorno=ww
   endelse     
   im[i,*,*]=tvrd()
   im0=im[i,*,*]                      
   device,/close                 
   set_plot,'x'        
   pp=where(reform(im[i,*,*]) eq 0)
   ppu=where(reform(im[i,*,*]) eq 255)
   if pp[0] ne -1 then im0[pp]=255 
   if ppu[0] ne -1 then im0[ppu]=0
   im[i,*,*]=im0    
   ;tv,w   
endfor                       
set_plot,'z'                  
DEVICE,set_RES=[4*250,2*250] 
for i=0,7 do begin
   x=(i mod 4)*250 & y=250-(i/4)*250
   tv,im[i,*,*],x,y
endfor     
fig=tvrd()
tvlct,r,g,b,/get
;if not keyword_set(ps) then device,/close
device,/close
;set_plot,'x'
;for i=0,nf-1 do begin
;   show_loop,nome=nome,ind=i
;endfor
rad='/Users/usuario/Dropbox/ProgramasIDL/Tratamento_NORH/Figuras/'
p=strpos(nome,'nn')
q=strpos(nome,'.sav')
nn=rad+'int'+strmid(nome,p+3,q-p-2)+'.png'
if keyword_set(ps) then begin
;   mydevice = !D.NAME
   ; Set plotting to PostScript:
   SET_PLOT, 'PS'
;   loadct,39,/silent
   tvlct,r,g,b
   tv,fig,0.41,0.095+DC,/norm,xsize=0.57,ysize=.4
endif else begin
   write_png,nn,fig,r,g,b
endelse
set_plot,mydevice
end
