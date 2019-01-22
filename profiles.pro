pro profiles,nome
restore,nome
;window,0,xsiz=1200,ysize=300
;!p.multi=[0,4,1]
mydevice=!d.name
set_plot, 'PS' 
pp=where(mag_res.b ne 0)
dx=0.05
dy=0.08



plot,mag_res.b[pp],mag_res.raio_v[pp],psy=2,back=255,col=0,/yno,chars=.5,xtit='Magnetic Field (G)',ytit='Height (Solar Radius)',position=[0.+dx,0+dy,0.15+dx,0.2+dy],/noer,syms=0.3
plot,mag_res.np[pp],mag_res.raio_v[pp],/yno,psy=2,back=255,col=0 ,chars=.5,xtit='Ambient Density (cm-3)', ytit='Height (Solar Radius)',position=[0.20+dx,0+dy,0.35+dx,0.2+dy],/noer,syms=0.3,/xlog
plot,20.*mag_res.np[pp]/mag_res.b[pp]/1.e9,mag_res.raio_v[pp],psy=2,/yno,/nodat,back=255,col=0,chars=.5,xtit='Razin Frequency (GHz)', ytit='Height (Solar Radius)',position=[0.4+dx,0+dy,0.55+dx,0.2+dy],/noer,syms=0.3,xticks=4
oplot,20.*mag_res.np[pp]/mag_res.b[pp]/1.e9,mag_res.raio_v[pp],psy=2,col=180,syms=0.3
plot,mag_res.nel[pp],mag_res.comp[pp],yran=[1.,0.],psy=2,back=255,col=0,chars=.5,xtit='Non-Thermal Density(cm-3)',ytit='Length from top (normalized)',position=[0.59+dx,0+dy,0.74+dx,0.2+dy],/noer,syms=0.3,xra=[1.e6,max(mag_res.nel[pp])],xticks=3,/xlog
plot,mag_res.theta[pp],mag_res.comp[pp],yran=[1.,0.],psy=2,back=255,col=0,chars=.5,xtit='View Angle(deg)',ytit='Length from top (normalized)',position=[0.78+dx,0+dy,0.93+dx,0.2+dy],/noer,syms=0.3,xra=[0.,180],xst=1,xminor=3,xticks=4

end
