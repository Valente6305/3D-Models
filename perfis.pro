pro perfis,nome
window,0,xsiz=1200,ysize=300
!p.multi=[0,4,1]
pp=where(mag_res.b ne 0)
plot,mag_res.b[pp],mag_res.raio_v[pp],psy=3,back=255,col=0,/yno,chars=2.,xtit='Magnetic Field (G)',ytit='Height (Solar Radius)'
plot,mag_res.np[pp],mag_res.raio_v[pp],/yno,psy=3,back=255,col=0 ,chars=2.,xtit='Ambient Density (cm-3)', ytit='Height (Solar Radius)'
plot,20.*mag_res.np[pp]/mag_res.b[pp]/1.e9,mag_res.raio_v[pp],psy=3,/yno,/nodat,back=255,col=0,chars=2.,xtit='Razin Frequency (GHz)', ytit='Height (Solar Radius)'
oplot,20.*mag_res.np[pp]/mag_res.b[pp]/1.e9,mag_res.raio_v[pp],psy=3,col=180
plot,mag_res.nel[pp],mag_res.comp[pp],yran=[0.,2],psy=3,back=255,col=0,chars=2.,xtit='Non-Thermal Density(cm-3)',ytit='Length from top (normalized)'

end
