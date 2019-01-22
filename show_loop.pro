function m_transf, st_arc,problem
  ; Teste do diametro do arco (para nao conter o eixo mag.)
  ; e definicao do b_arc e matriz de transformacao
  ; A .tr que retorna nesta routina serve para transformar 
  ; matrizes que estao no sistema de coord magnetico com 
  ; fase=0 (spot voltado para o obs) e com visao equatorial (mag.).
  ; O resultado apos a transformacao e' o spot voltado para 
  ; o obs com visao equatorial (equador de giro).
;
; problema nao existe de inicio
;
  problem=0
;  if (90.-st_arc.b_lat) gt 90. then beta_norte=180.-(90.0-st_arc.B_lat) $
;  else beta_norte=(90.0-st_arc.B_lat)
;  dd=abs(sphdist(0.0,st_arc.B_lat,0.0,90.-beta_norte, $
;            /degree))
;  st_arc.b_arc=dd

;  a=(90.-st_arc.b_lat)*!dtor 
;  b=st_arc.b_arc*!dtor 
;  c=(90.-st_arc.b_lat)*!dtor
;  rotat=acos(((cos(a)-cos(b)*cos(c))/(sin(b)*sin(c))<1.)>(-1.0))/!dtor
;  if a lt 0.0000001 then rotat=0.0
;  rotat=-rotat
;  T3D, /Reset, Rotate=[-st_arc.b_arc,rotat,0.],matrix=tmp
;  T3D, tmp, Rotate=[(90.0-st_arc.b_lat),0.,0.],matrix=tmp ; Tirei o
;  PHI
  T3D,/Reset, Rotate=[0.,0.,90.-st_arc.azim],matrix=tmp
  T3D, tmp, Translate=[0.,0.,(st_arc.l_arc-st_arc.deep)*st_arc.r_sun],matrix=tmp 
  T3D, tmp, Rotate=[-st_arc.B_lat,0.,0.],matrix=tmp 
  return, tmp
  erro:
end

pro show_loop,st_arc,px,ind=ind, graf=graf, nome=nome,contorno=contorno
 A = FINDGEN(17) * (!PI*2/16.)
 USERSYM, COS(A), SIN(A) ;, /FILL
;=================================================================
; plot,lim.coord[*,*,0]*960.,lim.coord[*,*,1]*960.,/nodata,/iso
; for i=0,63 do oplot,lim.coord[i,*,0]*960.,lim.coord[i,*,1]*960.
      if keyword_set(nome) then begin
        restore,nome,/v
        pp=strpos(nome,'phi')
        name=strmid(nome,pp,255)
        st_arc=mag_res.st_arc
        px=mag_res.px
      endif
      spx=st_arc.px

     if keyword_set(graf) then begin
         spx=0.0
         lpx=1.0
      endif else lpx=px
      v_angle=90.-st_arc.i 
      z=st_arc.l_arc*cos(st_arc.b_lat*!dtor)
      y=st_arc.l_arc*sin(st_arc.b_lat*!dtor)
      x=0.0
      T3D,/reset,rotate=[0.0,st_arc.phi,0.0],matrix=tmp
      T3D,tmp,rotate=[v_angle,0.0,0.0],matrix=tmp
      cen=[[x],[y],[z],[1.]]#tmp
      limite_arc,st_arc,lim,n_linh=64
      imagemc=intarr(37,300,400)
      if keyword_set(graf) then begin
         set_plot,'X'
         window, xs=300, ys=400,/free
         loadct,3,/silent
         tvlct,RRR,GGG,BBB,/get
         !p.position=[0,0,1,1]
         set_isoxy, cen[0,0]+lim.x_min*2,cen[0,0]-lim.x_min*2, $
                    cen[0,1]+lim.x_min*2,cen[0,1]-lim.x_min*2  
         limite_arc,st_arc,lim,n_linh=64
      endif
      T3D, /reset, Rotate=[v_angle,0.,0.],matrix=tmp
      T3D, /reset, Rotate=[st_arc.i-90.+v_angle,0.,0.],matrix=tmp2
      for w=0,n_elements(lim.coord[*,0,0])-1l do lim.coord[w,*,*]= $
                                            reform(lim.coord[w,*,*])#tmp
;=========================================
      lim.spot=lim.spot#tmp
;=========================================
      cc=[[lim.x0],[lim.y0],[lim.z0],[1]]#tmp2
      cc2=[[lim.dip_coord[0]],[lim.dip_coord[1]],[lim.dip_coord[2]],[1.0]]#tmp2 
      st_arc.tr=m_transf(st_arc,problem)
      if (problem eq 1) then begin
	  status=-1
	  goto, erro
       endif
      if keyword_set(graf) then begin
         sphinit, long=0, lat=v_angle, rad=1, fill=200, pa=0;  245, pa=0
         for lng=-180,180,10 do sphlng,lng+st_arc.phi,1
         sphlng,st_arc.phi,1, col=100, thic=2
         for lat=-80,80,10 do sphlat,lat,1
         sphlat,0.,1, col=100, thic=2
         sphinit, long=0, lat=v_angle, rad=1, pa=0
         sphrad,0.,90.,1.,1.55
      endif
      if keyword_set(file) then www=tvrd(/true)
      if keyword_set(file) then begin 
         www[0,*,*]=rotate(reform(www[0,*,*]),7)
         www[1,*,*]=rotate(reform(www[1,*,*]),7)
         www[2,*,*]=rotate(reform(www[2,*,*]),7)
      endif
;=======================================================================
; parte do grayplot
if keyword_set(nome) then begin
   if n_elements(nome) eq 0 then name='Spectrum'
   n=n_elements(imagem.image)
   sz=size(imagem.image[*].imagem[*,*,*])
   stokes_I=fltarr(n,sz[2],sz[3])
   stokes_V=stokes_I
   fluxos=fltarr(n)
   nsq=ceil(sqrt(n))
   spix=128
   fmt='('+strcompress(string(n)+'f5.1)',/remove_all)
;   tit=string(imagem.freq/1.e9,format=fmt)+' GHz'
   tit=strcompress(string(imagem.freq/1.e9)) 
;   window, xs=spix*nsq, ys=spix*nsq,/free,title=tit 
   for i=0, n-1 do begin 
      ii=i/nsq
      jj=i mod nsq
      stokes_I[i,*,*]=smooth(imagem.image[i].imagem[0,*,*]+imagem.image[i].imagem[1,*,*],1)
      stokes_V[i,*,*]=imagem.image[i].imagem[0,*,*]-imagem.image[i].imagem[1,*,*]
      im=bytscl(rebin(reform(stokes_I[i]),spix,spix))
   ;if i eq 0 then save,file='imagem.sav',im,spix,stokes_I,mag_res
 ;     tv,im,jj*spix,(nsq-ii-1)*spix 
 ;     fluxos[i]=total(stokes_I)
 ;     xyouts,jj*spix+fix(0.39*spix),(nsq-ii-1)*spix+fix(0.937*spix),fluxos[i],/dev
   endfor
endif
n=mag_res.n
delta=2.*abs(lim.x_min)/(n-1.) ;*mag_res.px/st_arc.px; arcsec em raio da estrela
x0=lim.center[0,0] ;x do centro do sol em raio da estrela
y0=lim.center[0,1] ;y do centro do sol
spc_min=-(n-1-(n-1)/2.)*delta
spc_max=(n-1-(n-1)/2.)*delta
arr=(findgen(n)-(n-1.)/2.)*delta
grayplot,reform(stokes_I[ind,*,*]),xarr=(arr+x0)*960., yarr=(arr+y0)*960.0,string=tit[ind]+' GHz' ; arcsec
;grayplottv,bytscl(reform(stokes_I)),xarr=(arr+x0)*960.,
;yarr=(arr+y0)*960.0 ; arcsec
if keyword_Set(contorno) then contour,contorno,(arr+x0)*960., (arr+y0)*960.0, xst=4,yst=4,levels=max(contorno)*[0.5],/over,col=255,thick=3.0
if keyword_Set(contorno) then contour,contorno,(arr+x0)*960., (arr+y0)*960.0, xst=4,yst=4,levels=max(contorno)*[0.5],/over,col=0
;================================fim do grayplot==============================
      for w=0,n_elements(lim.coord[*,0,0])-1l,4 do begin ; Apenas 8 arcos, OK!!
           dd=sqrt(lim.coord[w,*,0]^2+lim.coord[w,*,1]^2)
           pp=where(lim.coord[w,*,2] ge 0 or dd ge 1.,compl=npp)
           if npp[0] ne -1 then lim.coord[w,npp,*]=!values.f_nan
           oplot,(lim.coord[w,*,0]-0*(lim.center[0,0]-lim.x0))*960., $
                 (lim.coord[w,*,1]-0*(lim.center[0,1]-lim.y0))*960.,thick=1,col=80 ;100 
        endfor

      plots,cc[0,0],cc[0,1],psy=3,col=255 ;50
      plots,cc2[0,0],cc2[0,1],psy=8,col=255
;stop
erro:
      if keyword_set(graf) then begin
         set_isoxy,0,0,0,0
         !p.position=0
      endif
end
