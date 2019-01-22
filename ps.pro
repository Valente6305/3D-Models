pro ps,name,file=file
;INPUT:  name é o ff do file search
;        file é o nome do pdf (nao incluir a extensao)
;=================== Setup inicial=======================
mydevice = !D.NAME
; Set plotting to PostScript:
SET_PLOT, 'PS'

if keyword_set(file) then begin
   pp=strpos(file,'.')
   if pp[0] ne -1 then file=strmid(file,0,pp)
endif 

if keyword_set(file) then filename=file else $
   filename='~/postscript_geral' ; SEM O .PS !!!!!!

DC=0.22

; Use DEVICE to set some PostScript device options:
DEVICE, FILENAME=filename+'.ps', /LANDSCAPE,/color,xoffset=0.2,yoffset=25

;Leitura do fundo
back=read_png('Slide1.png')
nnn=n_elements(name)

for iii=0l, nnn-1l do begin
nome=name[iii]
;DEVICE, GET_DECOMPOSED=orig_decomposed, DECOMPOSED=0
LOADCT, 0,/silent

; Make a simple plot to the PostScript file:
tv,back,0.,0.+DC,/norm,xsize=1,ysize=1.,true=1
loadct,3,/silent
;================== End Setup Inicial====================

;################# Coloque os plots aqui ################

tvtotal,nome,/ps
show_sois,nome,/ps
show_loops,nome,/ps
profiles,nome
roda_imag,nome,/ps
p=strpos(nome,'nn')
q=strpos(nome,'.sav')
nn='nn'+strmid(nome,p+2,q-p-2)
xyouts,0.05,1.13,nn,/norm,chars=1.0,charthi=2

;PLOT, FINDGEN(10),position=[0.0,0.0,.3,.4],/noerase
;Plot,findgen(10)+30,positi=[0.7,0.0,0.99,0.99],/noerase
;tv,bytscl(dist(50)),0.0,0.5,/norm,xsiz=0.5,ysiz=0.3
;tvlct,r,g,b
;tv,back,0.,0.,/norm,xsize=1,ysize=1.,true=1
plot,findgen(2),/nodata,xst=4,yst=4
endfor

; ########## Fim da area de plotagens ##################
;==================== Setup Final ======================
; Close the PostScript file:
set_plot,'PS'
DEVICE, /CLOSE
; Return plotting to the original device:
SET_PLOT, mydevice
;================== End Setup Final ====================

spawn,'ps2pdf -dAutoRotatePages=/None '+filename+'.ps '+filename+'.pdf'
spawn,'/usr/local/bin/pdftk '+filename+'.pdf cat 1-endL output '+filename+'h.pdf'

end
