pro tvtotal,nome,ps=ps
restore,nome
;mag_res
dc=0.22
dx=50
dy=20
cor_back=250
cor_l=0
loadct,3,/silent
;if keyword_set(ps) then begin
;;   window,0,xsiz=3*64+dx,ysiz=4*64+dy
;   device,set_res=[3*64+dx,5*64+dy]
;endif else begin
   set_plot,'Z'
   DEVICE,set_RES=[3*64+dx,5*64+dy] 
;endelse     
tvlct,R,G,B,/get

tv,bytarr(3*64+dx,5*64+dy)+255

var=mag_res.b
ppp=where(var ne 0)
bmin=min(var[ppp],max=bmax)
bmed=mean(var[ppp])
x=[0,64,128]+dx
y=[128,128,128]+2*64+dy
im1=bytscl(rebin(total(var,1),64,64))
pp1=where(im1 eq 0)
contour,im1,findgen(64),findgen(64),levels=[0.99],path_xy=xy1,/path_data_coords
im1=bytscl(rebin(max(var,dimen=1),64,64))
im1[reform(xy1[0,*]),reform(xy1[1,*])]=cor_l
im1[0,*]=255
im1[63,*]=255
im1[*,0]=255
im1[*,63]=255
im1[pp1]=cor_back
tv,im1,x[0],y[0]
im2=bytscl(rebin(total(var,2),64,64))
pp2=where(im2 eq 0)
contour,im2,findgen(64),findgen(64),levels=[0.99],path_xy=xy2,/path_data_coords
im2=bytscl(rebin(max(var,dimen=2),64,64))
im2[reform(xy2[0,*]),reform(xy2[1,*])]=cor_l
im2[0,*]=255
im2[63,*]=255
im2[*,0]=255
im2[*,63]=255
im2[pp2]=cor_back
tv,im2,x[1],y[1]
im3=bytscl(rebin(total(var,3),64,64))
pp3=where(im3 eq 0)
contour,im3,findgen(64),findgen(64),levels=[0.99],path_xy=xy3,/path_data_coords
im3=bytscl(rebin(max(var,dimen=3),64,64))
im3[reform(xy3[0,*]),reform(xy3[1,*])]=cor_l
im3[0,*]=255
im3[63,*]=255
im3[*,0]=255
im3[*,63]=255
im3[pp3]=cor_back
tv,im3,x[2],y[2]

var=mag_res.nel
nelmin=min(var[ppp],max=nelmax)
nelmed=mean(var[ppp])
x=[0,64,128]+dx
y=[64,64,64]+2*64+dy
;im1=bytscl(rebin(total(var,1),64,64))
im1=bytscl(rebin(max(var,dimen=1),64,64))
im1[reform(xy1[0,*]),reform(xy1[1,*])]=cor_l
im1[0,*]=255
im1[63,*]=255
im1[*,0]=255
im1[*,63]=255
im1[pp1]=cor_back
tv,im1,x[0],y[0]
;im2=bytscl(rebin(total(var,2),64,64))
im2=bytscl(rebin(max(var,dimen=2),64,64))
im2[reform(xy2[0,*]),reform(xy2[1,*])]=cor_l
im2[0,*]=255
im2[63,*]=255
im2[*,0]=255
im2[*,63]=255
im2[pp2]=cor_back
tv,im2,x[1],y[1]
;im3=bytscl(rebin(total(var,3),64,64))
im3=bytscl(rebin(max(var,dimen=3),64,64))
im3[reform(xy3[0,*]),reform(xy3[1,*])]=cor_l
im3[0,*]=255
im3[63,*]=255
im3[*,0]=255
im3[*,63]=255
im3[pp3]=cor_back
tv,im3,x[2],y[2]

var=mag_res.np
npmin=min(var[ppp],max=npmax)
npmed=mean(var[ppp])
x=[0,64,128]+dx
y=[0,0,0]+2*64+dy
;im1=bytscl(rebin(total(var,1),64,64))
im1=bytscl(rebin(max(var,dimen=1),64,64))
im1[reform(xy1[0,*]),reform(xy1[1,*])]=cor_l
im1[0,*]=255
im1[63,*]=255
im1[*,0]=255
im1[*,63]=255
im1[pp1]=cor_back
tv,im1,x[0],y[0]
;im2=bytscl(rebin(total(var,2),64,64))
im2=bytscl(rebin(max(var,dimen=2),64,64))
im2[reform(xy2[0,*]),reform(xy2[1,*])]=cor_l
im2[0,*]=255
im2[63,*]=255
im2[*,0]=255
im2[*,63]=255
im2[pp2]=cor_back
tv,im2,x[1],y[1]
;im3=bytscl(rebin(total(var,3),64,64))
im3=bytscl(rebin(max(var,dimen=3),64,64))
im3[reform(xy3[0,*]),reform(xy3[1,*])]=cor_l
im3[0,*]=255
im3[63,*]=255
im3[*,0]=255
im3[*,63]=255
im3[pp3]=cor_back
tv,im3,x[2],y[2]



var=mag_res.T
tmin=min(var[ppp],max=tmax)
tmed=mean(var[ppp])
x=[0,64,128]+dx
y=[0,0,0]+64+dy
;im1=bytscl(rebin(total(var,1),64,64))
im1=bytscl(rebin(max(var,dimen=1),64,64))
im1[reform(xy1[0,*]),reform(xy1[1,*])]=cor_l
im1[0,*]=255
im1[63,*]=255
im1[*,0]=255
im1[*,63]=255
im1[pp1]=cor_back
tv,im1,x[0],y[0]
;im2=bytscl(rebin(total(var,2),64,64))
im2=bytscl(rebin(max(var,dimen=2),64,64))
im2[reform(xy2[0,*]),reform(xy2[1,*])]=cor_l
im2[0,*]=255
im2[63,*]=255
im2[*,0]=255
im2[*,63]=255
im2[pp2]=cor_back
tv,im2,x[1],y[1]
;im3=bytscl(rebin(total(var,3),64,64))
im3=bytscl(rebin(max(var,dimen=3),64,64))
im3[reform(xy3[0,*]),reform(xy3[1,*])]=cor_l
im3[0,*]=255
im3[63,*]=255
im3[*,0]=255
im3[*,63]=255
im3[pp3]=cor_back
tv,im3,x[2],y[2]


var=mag_res.theta
thetamin=min(var[ppp],max=thetamax)
thetamed=mean(var[ppp])
cc=where(var gt 90)
if cc[0] ne -1 then var[cc]=180.0-var[cc]
x=[0,64,128]+dx
y=[0,0,0]+dy
;im1=bytscl(rebin(total(var,1),64,64))
im1=bytscl(rebin(max(var,dimension=1),64,64))
im1[reform(xy1[0,*]),reform(xy1[1,*])]=cor_l
im1[0,*]=255
im1[63,*]=255
im1[*,0]=255
im1[*,63]=255
im1[pp1]=cor_back
tv,im1,x[0],y[0]
;im2=bytscl(rebin(total(var,2),64,64))
im2=bytscl(rebin(max(var,dimen=2),64,64))
im2[reform(xy2[0,*]),reform(xy2[1,*])]=cor_l
im2[0,*]=255
im2[63,*]=255
im2[*,0]=255
im2[*,63]=255
im2[pp2]=cor_back
tv,im2,x[1],y[1]
;im3=bytscl(rebin(total(var,3),64,64))
im3=bytscl(rebin(max(var,dimen=3),64,64))
im3[reform(xy3[0,*]),reform(xy3[1,*])]=cor_l
im3[0,*]=255
im3[63,*]=255
im3[*,0]=255
im3[*,63]=255
im3[pp3]=cor_back
tv,im3,x[2],y[2]

xyouts,3,4*64+32+dy,'B!Imax!N',size=1.4,charth=1.7,/device,col=0
xyouts,3,3*64+32+dy,'Nel!Imax!N',size=1.4,charth=1.7,/device,col=0
xyouts,3,2*64+32+dy,'Np!Imax!N',size=1.4,charth=1.7,/device,col=0
xyouts,3,1*64+32+dy,'T!Imax!N',size=1.4,charth=1.7,/device,col=0
xyouts,3,0*64+32+dy,'View!Imax!N',size=1.4,charth=1.7,/device,col=0

xyouts,0*64+32+dx,1,'X',size=1.5,charth=2.0,/device,col=0
xyouts,1*64+32+dx,1,'Y',size=1.5,charth=2.0,/device,col=0
xyouts,2*64+32+dx,1,'Z',size=1.5,charth=2.0,/device,col=0


im=tvrd()
rad='/Users/usuario/Dropbox/ProgramasIDL/Tratamento_NORH/Figuras/tot/'
p=strpos(nome,'nn')  ; era img
q=strpos(nome,'.sav')
;fig=tvrd()
nn=rad+'tot'+strmid(nome,p+2,q-p-2)+'.png'    ; era p+3 (devido ao img)
if not keyword_set(ps) then write_png,nn,im,r,g,b
fmt='(2E9.2)'
if keyword_set(ps) then begin
   set_plot,'PS'
   dc0=0.03 ; MIN
   Emin=imagem.gyro.energy[0]
   Emax=imagem.gyro.energy[1]
   delta=imagem.gyro.delta
   Emed=(Emax^(-delta+2)-(Emin^(-delta+2)))/(Emax^(-delta+1)-(Emin^(-delta+1)))*(1.-delta)/(2.-delta)
   tv,im,0.198,0.097+dc,/norm,xsize=0.199,ysize=.33
   xyouts,0.58+dc0,0.9800,string(bmin,format='(f8.0)'),/norm,chars=1.2,charthi=3
   xyouts,0.58+dc0,0.9515,string(nelmin,format=fmt),/norm,chars=1.2,charthi=3
   xyouts,0.58+dc0,0.9230,string(npmin,format=fmt),/norm,chars=1.2,charthi=3
   xyouts,0.58+dc0,0.8945,string(tmin,format=fmt),/norm,chars=1.2,charthi=3
   xyouts,0.58+dc0,0.8660,string(thetamin,format='(f6.0)'),/norm,chars=1.2,charthi=3
   xyouts,0.58+dc0,0.8375,string(Emin,format='(f6.0)'),/norm,chars=1.2,charthi=3
   xyouts,0.60+dc0,0.8090,string(delta,format='(f4.1)'),/norm,chars=1.2,charthi=3

   dc=0.12 ; MAX
   xyouts,0.58+dc0+dc,0.9800,string(bmax,format='(f8.0)'),/norm,chars=1.2,charthi=3
   xyouts,0.58+dc0+dc,0.9515,string(nelmax,format=fmt),/norm,chars=1.2,charthi=3
   xyouts,0.58+dc0+dc,0.9230,string(npmax,format=fmt),/norm,chars=1.2,charthi=3
   xyouts,0.58+dc0+dc,0.8945,string(Tmax,format=fmt),/norm,chars=1.2,charthi=3
   xyouts,0.58+dc0+dc,0.8660,string(thetamax,format='(f6.0)'),/norm,chars=1.2,charthi=3
   xyouts,0.58+dc0+dc,0.8375,string(Emax,format='(f7.0)'),/norm,chars=1.2,charthi=3
   xyouts,0.60+dc0+dc,0.8090,string(delta,format='(f4.1)'),/norm,chars=1.2,charthi=3

   dc=0.12*2 ; MED
   xyouts,0.58+dc0+dc,0.9800,string(bmed,format='(f8.0)'),/norm,chars=1.2,charthi=3
   xyouts,0.58+dc0+dc,0.9515,string(nelmed,format=fmt),/norm,chars=1.2,charthi=3
   xyouts,0.58+dc0+dc,0.9230,string(npmed,format=fmt),/norm,chars=1.2,charthi=3
   xyouts,0.58+dc0+dc,0.8945,string(tmed,format=fmt),/norm,chars=1.2,charthi=3
   xyouts,0.58+dc0+dc,0.8660,string(thetamed,format='(f6.0)'),/norm,chars=1.2,charthi=3
   xyouts,0.58+dc0+dc,0.8375,string(Emed,format='(f7.1)'),/norm,chars=1.2,charthi=3
;   xyouts,0.58+dc,0.8090,string(bmin),/norm,chars=1.2,charthi=3

; Parte da tabela 1
   dy=0.06
   pqp=strpos(ofile,'ie')
   ww=fix(strmid(ofile,pqp+2,4))
   spar='Homog'
   if ww eq 5 then spar='Feet '
   if ww eq 200 then spar='Top  ' 

   xyouts,0.07+dc0,0.4800+dy,string(mag_res.st_arc.b_lat,format='(f5.1)'),/norm,chars=1.2,charthi=3
   xyouts,0.07+dc0,0.4515+dy,string(mag_res.st_arc.phi,format='(f5.1)'),/norm,chars=1.2,charthi=3
   xyouts,0.07+dc0,0.4230+dy,string(mag_res.st_arc.azim,format='(f5.1)'),/norm,chars=1.2,charthi=3
   xyouts,0.07+dc0,0.3945+dy,string(mag_res.st_arc.h_arc,format='(f6.3)'),/norm,chars=1.2,charthi=3
   xyouts,0.07+dc0,0.3660+dy,string(mag_res.st_arc.r_arc,format='(f6.3)'),/norm,chars=1.2,charthi=3
   xyouts,0.07+dc0,0.3375+dy,string(mag_res.st_arc.foot_s,format='(f6.3)'),/norm,chars=1.2,charthi=3
   xyouts,0.08+dc0,0.3090+dy,spar,/norm,chars=1.2,charthi=3



endif
end
