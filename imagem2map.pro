PRO imagem2map,imagem,map_i,map_v,px2cm=px2cm

 ;; routine to create SSW map structures using output from calc_image
 ;; Written: 02-May-2012, PJSimoes, Glasgow
 ;; NOTE: if no date/time is set, set B0=0.0 when using plot_map
  
 arc2cm=72527093.
 cm2arc=1./arc2cm

 data=fltarr(imagem.n,imagem.n)
 latlon=[imagem.st_arc.b_lat,imagem.st_arc.phi]
 px2cm=mag_res.px 
 px=px2cm*cm2arc
 freq=imagem.freq

 ;; SSW routines
 pos = hel2arcmin(latlon[0],latlon[1],B0=0.0,date=date) * 60.0
 ;;---------------

 xc=pos[0]
 yc=pos[1]
 dx=px
 dy=px
 xunits='arcsec'
 yunits='arcsec'
 dur=0.
 ID = 'gyro3D '+string(freq/1e9,format='(f0.1)')+' GHz'
 roll_center = 0.0
 roll_angle = 0.0

 ;; map struct
 map={data:data,xc:xc,yc:yc,dx:dx,dy:dy,time:'',id:'' $
      ,dur:dur,xunits:xunits,yunits:yunits,roll_angle:roll_angle $
      ,roll_center:roll_center}

 map_i=replicate(map,n_elements(freq))
 map_v=map_i

 FOR i=0,n_elements(freq)-1 DO BEGIN  
    a=reform(imagem.image[i].imagem[0,*,*]+imagem.image[i].imagem[1,*,*])
    b=reform(imagem.image[i].imagem[0,*,*]-imagem.image[i].imagem[1,*,*])
    map_i[i].data=a
    map_V[i].data=b
    map_i[i].id=id[i]+' I'
    map_v[i].id=id[i]+' V'
 ENDFOR 
END 
;;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PRO test_helio_pos,lat,lon

 loadct,0
 lat=imagem.st_arc.b_lat
 lon=imagem.st_arc.phi
 b0=0.0
 n=32
 data=bytarr(n,n)
 
 x=gaussian(findgen(n),[1.0,n/2.,n/8.])
 y=gaussian(findgen(n),[1.0,n/2.,n/8.])
 data=x#y

 pnt=hel2arcmin(lat,lon,b0=b0)*60.0
 map=make_map(data,b0=b0,xc=pnt[0],yc=pnt[1],dx=5,dy=5)
!p.multi=[0,2,1]
 plot_map,map,mark=pnt $
          ,center=[0,0],fov=40,grid=10 $
          ,b0=0,gsty=0,gcol=150
 plot_map,map,grid=5,fov=5
!p.multi=0

END 
