;-------------------------------------------------------------
;+
; NAME:
;       SPHDIST
; PURPOSE:
;       Angular distance between points on a sphere.
; CATEGORY:
; CALLING SEQUENCE:
;       d = sphdist(long1, lat1, long2, lat2)
; INPUTS:
;       long1 = longitude of point 1.         in
;       lat1 = latitude of point 1.           in
;       long2 = longitude of point 2.         in
;       lat2 = latitude of point 2.           in
; KEYWORD PARAMETERS:
;       Keywords:
;         /DEGREES means angles are in degrees, else radians.
; OUTPUTS:
;       d = angular distance between points.  out
; COMMON BLOCKS:
; NOTES:
;       Notes: points 1 and 2 may be arrays.
; MODIFICATION HISTORY:
;       R. Sterner, 5 Feb, 1991
;       R. Sterner, 26 Feb, 1991 --- Renamed from sphere_dist.pro
;
; Copyright (C) 1991, Johns Hopkins University/Applied Physics Laboratory
; This software may be used, copied, or redistributed as long as it is not
; sold and this copyright notice is reproduced on each copy made.  This
; routine is provided as is without any express or implied warranties
; whatsoever.  Other limitations apply as described in the file disclaimer.txt.
;-
;-------------------------------------------------------------

        function sphdist, long1, lat1, long2, lat2, $
          help=hlp, degrees=degrees

        if (n_params(0) lt 4) or keyword_set(hlp) then begin
          print,' Angular distance between points on a sphere.'
          print,' d = sphdist(long1, lat1, long2, lat2)'
          print,'   long1 = longitude of point 1.         in'
          print,'   lat1 = latitude of point 1.           in'
          print,'   long2 = longitude of point 2.         in'
          print,'   lat2 = latitude of point 2.           in'
          print,'   d = angular distance between points.  out'
          print,' Keywords:'
          print,'   /DEGREES means angles are in degrees, else radians.'
          print,' Notes: points 1 and 2 may be arrays.'
          return, -1
        endif

        cf = 1.0
        if keyword_set(degrees) then cf = !radeg

        ;--- Convert both points to rectangular coordinates. ---
        polrec, 1.0, lat1/cf, rxy, z1
        polrec, rxy, long1/cf, x1, y1
        polrec, 1.0, lat2/cf, rxy, z2
        polrec, rxy, long2/cf, x2, y2

        ;--- Compute vector dot product for both points. ---
        cs = x1*x2 + y1*y2 + z1*z2

        ;--- Compute the vector cross product for both points. ---
        xc = y1*z2 - z1*y2
        yc = z1*x2 - x1*z2
        zc = x1*y2 - y1*x2
        sn = sqrt(xc*xc + yc*yc + zc*zc)

        ;--- Convert to polar.  ------
        recpol, cs, sn, r, a
        return, cf*a

    end

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

  T3D,/Reset,Translate=[0.,0.,-st_arc.deep*st_arc.r_sun],matrix=tmp
  T3D, tmp, Rotate=[st_arc.d_ang,st_arc.a_inc,90.-st_arc.azim],matrix=tmp
  T3D, tmp, Translate=[0.,0.,st_arc.l_arc*st_arc.r_sun],matrix=tmp 
;  T3D, tmp, Translate=[0.,0.,st_arc.l_arc*st_arc.r_sun],matrix=tmp ; pto rotacao=superficie
  T3D, tmp, Rotate=[-st_arc.B_lat,0.,0.],matrix=tmp 

  return, tmp
  erro:
end

pro st_sun,st_arc,new=new,show=show,i=i,beta=beta,phi=phi, azim=azim,$
            mag=mag, dens=dens, px=px, r_sun=r_sun, nel=nel, $
            l_arc=l_arc, R_arc=R_arc,h_arc=h_arc, zoom=zoom, $
            x0=x0, y0=y0, b_lat=b_lat,print=print, $
            v_angle=v_angle, gira=gira, status=status,file=file, $
            foot_s=foot_s,a_inc=a_inc,assim=assim

; MODIFICATION HISTORY:
;       Criado em 2006 J.E.R.Costa
;       modificado em Outubro 2007 (eixo magnetico qqer longitude)
;       modificado junh0 2014 para incluir a_inc
;       modificado maio 2015 para incluir  d_ang (assimetria)
; **
; ** iniciando com variavel de problema zerada
;
status=0
time0=systime(/sec)  
;=================================================================
;                      Parametro Default
;=================================================================
if keyword_set(new) then begin
 st_arc={                  $ ; Definicao da Estrutura do Arco
   X0      :200.0        , $ ; Centro do Sol em X (pixel)
   Y0      :0.0          , $ ; Centro do Sol em Y (pixel)
   R_SUN   :100.0        , $ ; Raio do Sol em pixel (era 40)
   L_ARC   :1.002        , $ ; Altura da base do arco (em R_SUN) 1.05
   H_ARC   :0.02         , $ ; Altura do arco (em R_SUN) 0.1
   R_ARC   :0.01         , $ ; Raio do arco no topo (em R_SUN) 0.02
   Deep    :9.349e-3     , $ ; Profundidade do dipolo (em R_SUN) 0.047
;                              medido a partir da base L_ARC
;       Se h_arc=Foot_s => deep=[1/(1-cos^3(46.925))-1]*h_arc                   
   Foot_s  :0.02         , $ ; Separacao dos pes em R_SUN 0.1
   B_Lat   :30.          , $ ; Latitude do dipolo magnetico 
;   B_Long  :0.           , $ ; Longitude do dipolo magnetico
   Azim    :90.0         , $ ; Azimuth do dipolo em rel. ao equador
   PX      :6.960e08     , $ ; Escala cm/pixel (Sun_r(cm)/100)
   I       :90.0000      , $ ; Ang. do Obs. em rel. ao eixo rot.
   PHI     :30.0000      , $ ; Ang. de fase. do meridiano central.
   B       :1.00000e+03  , $ ; Campo magnetico na base do arco
   Nel     :1.00000e+7   , $ ; Densidade de eletrons nao termicos
   Np      :1.00000e+12  , $ ; Densidade de eletrons termicos
   B_arc   :90.0000      , $ ; Ang. entre o eixo mag e spot
   a_inc   :0.0          , $ ; Ang. de inclinacao do arco em graus
   d_ang   :0.0          , $ ; Ang. de assimetria dos pes do arco 
   assim   :0.0          , $ ; Fracao do aumneto de B=>(1.0-Bmaior/B)
   mirror  :15.8964      , $ ; Mirror ratio Bpes/Btopo 
   b0      :1590.31      , $ ; Constante magnetica (Gauss)
   tr      :fltarr(4,4)    $ ; Matriz de transformacao. Mag->Giro
   }
endif 
;=================================================================
 A = FINDGEN(17) * (!PI*2/16.)
 USERSYM, COS(A), SIN(A) ;, /FILL
;=================================================================
;                      Troca de Parametros
;=================================================================
   flag_deep=0
   if size(foot_s,/type) ne 0 or size(L_ARC,/type) ne 0 or $
      size(H_ARC,/type) ne 0 then begin
      if size(foot_s,/type) ne 0 then s=foot_s else s=st_arc.foot_s
      if size(L_arc,/type) ne 0 then L=l_arc else L=st_arc.L_arc
      if size(H_arc,/type) ne 0 then H=H_arc else H=st_arc.H_arc
      deep=findgen(351)/150.
      v=0.5*Sqrt(4*deep^2+s^2)-4*deep^2*(deep+H)/(4.*deep^2+s^2)
      if (where(v*v[1:*] le 0))[0] eq -1 then begin 
         print,string(7B),'Dipolo acima do pe ou linha muito tangente.'
         print,'Assumindo parametros default!'
         return
      endif
      flag_deep=1
      tmp=min(abs(v),ii)
      st_arc.deep=deep[ii]
      st_arc.L_arc=L
      st_arc.H_arc=H
      st_arc.Foot_s=s
;      st_arc.mirror=sqrt((st_arc.deep^2+st_arc.foot_s^2)/(4*st_arc.deep^2+$
 ;                st_arc.foot_s^2))*(4+(st_arc.foot_s/st_arc.deep)^2)^3/32.
   endif    
   if size(assim,/type) ne 0 then st_arc.assim=assim

   flag_rsun=0
   if size(R_SUN,/type) ne 0 then begin 
      st_arc.R_SUN=R_SUN
      flag_rsun=1
   endif
   flag_b=0
   if size(mag,/type) ne 0 then begin 
      st_arc.B=mag
      flag_b=1
   endif
   flag_b0=0
   if flag_rsun or flag_deep or flag_b then begin
     realh=st_arc.h_arc-st_arc.r_arc/2.
;     st_arc.B0=32.*st_arc.r_sun^3*st_arc.deep^6*sqrt(4*st_arc.deep^2+$
;        st_arc.foot_s^2)*$
;        (st_arc.deep+st_arc.h_arc)^3/(sqrt(st_arc.deep^2+st_arc.foot_s^2)*$
;        (4*st_arc.deep^2+st_arc.foot_s^2)^(3))*st_arc.b
     st_arc.B0=32.*st_arc.r_sun^3*st_arc.deep^6*sqrt(4*st_arc.deep^2+$
        st_arc.foot_s^2)*$
        (st_arc.deep+realh)^3/(sqrt(st_arc.deep^2+st_arc.foot_s^2)*$
        (4*st_arc.deep^2+st_arc.foot_s^2)^(3))*st_arc.b

r0=(st_arc.h_arc+st_arc.deep)*st_arc.r_sun
angf=atan(st_arc.foot_s/2,st_arc.deep) ; ang. polar da fotosfera
rf=r0*(cos(angf))^2                    ; raio vetor da fotosfera
print,r0,angf,st_arc.b0*sqrt(1.+3.*sin(angf)^2)/rf^3

     flag_b0=1
   endif
   if size(R_ARC,/type) ne 0 then begin 
      if R_ARC ge st_arc.h_arc then begin
      print,string(7B),'Seccao transversal muito gde!'
      st_arc.R_ARC=st_sun.h_arc
      endif else st_arc.R_ARC=R_ARC 
   endif  
   if size(i,/type) ne 0 then st_arc.i=i
   if size(phi,/type) ne 0 then st_arc.phi=phi
   if size(Azim,/type) ne 0 then st_arc.azim=azim
   if size(dens,/type) ne 0 then begin 
      st_arc.Np=dens
   endif
   if size(nel,/type) ne 0 then begin ; acrecentei nel e ambient na st_arc
      st_arc.nel=nel
   endif
   if size(px,/type) ne 0 then st_arc.px=px
   if size(R_SUN,/type) ne 0 then st_arc.R_SUN=R_SUN
      
   if size(x0,/type) ne 0 then st_arc.x0=x0 
   if size(y0,/type) ne 0 then st_arc.y0=y0 
   if size(B_LAT,/type) ne 0 then st_arc.B_LAT=B_LAT
   if size(a_inc,/type) ne 0 then st_arc.a_inc=a_inc
   if size(assim,/type) ne 0 or flag_deep or flag_rsun or flag_b0 then begin
      if size(assim,/type) eq 0 then assim=st_arc.assim 
      if assim lt 0 then signal_assim=-1.0 else signal_assim=1.0
      ;st_arc.assim=abs(assim)  eu tirei pq e' preciso saber
;      B0=32.*st_arc.r_sun^3*st_arc.deep^6*sqrt(4*st_arc.deep^2+st_arc.foot_s^2)*$
;      (st_arc.deep+st_arc.h_arc)^3/(sqrt(st_arc.deep^2+st_arc.foot_s^2)*$
;      (4*st_arc.deep^2+st_arc.foot_s^2)^(3))*st_arc.b

      r0=(st_arc.h_arc+st_arc.deep)*st_arc.r_sun
      angf=atan(st_arc.foot_s/2,st_arc.deep) ; ang. polar da fotosfera
      rf=r0*(cos(angf))^2                    ; raio vetor da fotosfera

;      k=((1./(1.+assim))*rf^3/sqrt(1+3.*sin(angf)^2))^2
;      r2=(findgen(301)/300.+0.5)*rf
;      v=r2^6-k*(1.+3.d0*(1.-r2/r0))
;      if (where(v*v[1:*] le 0))[0] eq -1 then stop
;      tmp=min(abs(v),ii)
;      r2=r2[ii]
;      ang2=acos(sqrt(r2/r0))
;      Bn=st_arc.b0*sqrt(1.0+3.*sin(ang2)^2)/r2^3
;      a2=r2^2+(st_arc.deep*st_arc.r_sun)^2-$
;         2.*r2*st_arc.deep*st_arc.r_sun*cos(ang2)
;      phi=!pi/2.-asin(r2*sin(ang2)/sqrt(a2))
;      st_arc.d_ang=signal_assim*phi/!dtor

      r2of=1.d0/(1.d0+abs(st_arc.assim))^(1./3)
      r2=r2of*rf
      ang2=acos(sqrt(r2of*cos(angf)^2))
      ctec=sqrt((st_arc.deep*st_arc.r_sun)^2+r2^2-2.*st_arc.deep*st_arc.r_sun*r2*cos(ang2))
      phitt=asin(r2*sin(ang2)/ctec)
      st_arc.d_ang=signal_assim*(!pi/2-phitt)/!dtor
  ; print,'era este d_ang=',st_arc.d_ang,assim
  ; stop
      if assim eq 0 then st_arc.d_ang=0.0 
   endif     

;   if size(B_LONG,/type) ne 0 then st_arc.B_LONG=B_LONG
   st_arc.tr=m_transf(st_arc,problem)
   if (problem eq 1) then begin
	status=-1
	goto, erro
   endif 

;O B_ang nao esta sendo usado. Deixei para manter o historico.
;   r_b=st_arc.r_sun*st_arc.b_arc
;   r_t=st_arc.r_sun*(st_arc.b_arc+st_arc.h_arc)
;   Da=asin(st_arc.r_arc/st_arc.b_arc)
   if size(v_angle,/type) eq 0 then v_angle=90.0-st_arc.i
;=================================================================

;=================================================================
;                         Parte Grafica
;=================================================================
   if keyword_set(show) then begin 
      z=st_arc.l_arc*cos(st_arc.b_lat*!dtor)
      y=st_arc.l_arc*sin(st_arc.b_lat*!dtor)
      x=0.0
      T3D,/reset,rotate=[0.0,st_arc.phi,0.0],matrix=tmp
      T3D,tmp,rotate=[v_angle,0.0,0.0],matrix=tmp
      cen=[[x],[y],[z],[1.]]#tmp
      limite_arc,st_arc,lim,n_linh=64
      imagem=intarr(37,300,400)
      set_plot,'X'
      window, xs=300, ys=400,/free
      loadct,3,/silent
      tvlct,RRR,GGG,BBB,/get
      !p.position=[0,0,1,1]
      if keyword_set(zoom) then begin
;         print,cen[0,0]+lim.x_min,cen[0,0]-lim.x_min, $
;                    cen[0,1]+lim.x_min,cen[0,1]-lim.x_min, $
;                    cen[0,0],lim.x_min, $
;                    cen[0,1]
         set_isoxy, cen[0,0]+lim.x_min*2,cen[0,0]-lim.x_min*2, $
                    cen[0,1]+lim.x_min*2,cen[0,1]-lim.x_min*2  
      endif else set_isoxy, -1.1, 1.1, -1.1, 1.1
      m_phi=st_arc.phi
      tphi0=0
      tphi9=0
      if keyword_set(gira) then begin
         tphi0=0
         tphi9=36
      endif
      for tphi=tphi0,tphi9 do begin
        erase, 10
        st_arc.phi=((tphi-tphi0)*10. + m_phi) mod 360.
        limite_arc,st_arc,lim,n_linh=64
        T3D, /reset, Rotate=[v_angle,0.,0.],matrix=tmp
        T3D, /reset, Rotate=[st_arc.i-90.+v_angle,0.,0.],matrix=tmp2
        for w=0,n_elements(lim.coord[*,0,0])-1l do lim.coord[w,*,*]= $
                                              reform(lim.coord[w,*,*])#tmp
;=========================================
        lim.spot=lim.spot#tmp
;stop
;=========================================
        cc=[[lim.x0],[lim.y0],[lim.z0],[1]]#tmp2
        cc2=[[lim.dip_coord[0]],[lim.dip_coord[1]],[lim.dip_coord[2]],[1.0]]#tmp2 
        st_arc.tr=m_transf(st_arc,problem)
        if (problem eq 1) then begin
	  status=-1
	  goto, erro
        endif
        sphinit, long=0, lat=v_angle, rad=1, fill=200, pa=0;  245, pa=0
        for lng=-180,180,10 do sphlng,lng+st_arc.phi,1
        sphlng,st_arc.phi,1, col=100, thic=2
        for lat=-80,80,10 do sphlat,lat,1
        sphlat,0.,1, col=100, thic=2
        sphinit, long=0, lat=v_angle, rad=1, pa=0
        sphrad,0.,90.,1.,1.55
;        sphrad,0.0+st_arc.phi,st_arc.b_lat,1., $
;               1.5,col=100,thic=3
        if keyword_set(file) then www=tvrd(/true)
        if keyword_set(file) then begin 
           www[0,*,*]=rotate(reform(www[0,*,*]),7)
           www[1,*,*]=rotate(reform(www[1,*,*]),7)
           www[2,*,*]=rotate(reform(www[2,*,*]),7)
        endif
        for w=0,n_elements(lim.coord[*,0,0])-1l do begin 
           dd=sqrt(lim.coord[w,*,0]^2+lim.coord[w,*,1]^2)
           pp=where(lim.coord[w,*,2] ge 0 or dd ge 1.,compl=npp)
           if npp[0] ne -1 then lim.coord[w,npp,*]=!values.f_nan
           oplot,lim.coord[w,*,0],lim.coord[w,*,1],thick=2,col=100 
        endfor
        plots,cc[0,0],cc[0,1],psy=3,col=255 ;50
        plots,cc2[0,0],cc2[0,1],psy=8,col=255
        wait,0.05
; ===============  Fim da linha do campo =============== 
         if keyword_set(file) then write_gif,file,tvrd(),RRR,GGG,BBB,/multiple
      endfor  
      st_arc.phi=m_phi
      set_isoxy,0,0,0,0
      !p.position=0
   endif      
;=================================================================

     if keyword_set(file) then write_gif,file,tvrd(),RRR,GGG,BBB,/close

; O return abaixo tem que ser o ultimo estagio do programa
     if keyword_set(print) then begin 
        ;help,st_arc;,/str
        return
     endif
  erro:
timef=systime(/sec)-time0 
print,'tempo de st_arc: ', timef 
end
