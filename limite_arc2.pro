pro limite_arc,st_sun,limite,n_phase=n_phase, n_linhas=n_linhas
 if size(n_phase,/type) eq 0 then n_phase=16
 if size(n_linhas,/type) eq 0 then n_l=32 else n_l=n_linhas
 ; Trinta e duas linhas no entorno da elipse.
 ; 127 pontos por linha.
 ang=(reverse((findgen(127)-63.)/63.)) ; -1< ang < +1
 coordz=fltarr(n_l+1,127)
 coordy=fltarr(n_l+1,127)
 coordx=fltarr(n_l+1,127)              ; todos zeros
 coord1=fltarr(n_l+1,127)+1.
 coord=fltarr(n_phase,n_l+1,127,4)

 r_ang=findgen(n_l)/float(n_l-1)*360.*!dtor
; burst_ang=atan(st_sun.r_arc,st_sun.l_arc+st_sun.h_arc-st_sun.deep)/!dtor
 burst_ang=atan(st_sun.r_arc,st_sun.deep+st_sun.h_arc)/!dtor
 zz=burst_ang*cos(r_ang) ; offsets de long. em graus
 yy=sqrt(1.0-(zz/burst_ang)^2)*st_sun.R_Arc ; ok
 yy[n_l/2:*]=-yy[n_l/2:*]
 zz=[zz,0.0]
 yy=[yy,0.0]+st_sun.deep+st_sun.H_arc
 r0=yy                ; ctes das linhas da seccao circular
 lamb=findgen(301)/300.*90.*!dtor
 L=st_sun.L_arc
 D=st_sun.L_arc-st_sun.deep  ; distancia do centro do sol ao dipolo
 for i=0,n_l do begin
    v=-D*r0[i]*cos(lamb)*cos(2.*lamb-!dpi)+r0[i]^2*(cos(lamb))^4+ $
       D*r0[i]*cos(lamb)-L^2+D^2
    if (where(v*v[1:*] le 0))[0] eq -1 then begin 
       print,string(7B),'Problema do calculo de lat. mag.'
       stop
   endif
   tmp=min(abs(v),ii)
   yy[i]=lamb[ii]
 endfor
 ; ====== Calculo das n_l linhas que saem da elipse da base ===========
 ; Aqui tem se definir a nova base do loop apos a assimetria d_ang !!!
for i=0,n_l-1+1 do begin
    d_ang=ang*yy[i]
    coordz[i,*]=r0[i]*((cos(d_ang))^3)
    coordy[i,*]=r0[i]*((cos(d_ang))^2)*sin(d_ang) 
 endfor
; =====================================================================
; =====================================================================
                           ; Coord esta' no sistema de
                           ; coordenadas magnetico planado no
                           ; limbo meridiano central. Visao Equatorial.
                           ; ==========================================
 m_lamb=atan(coordy[n_l,63],coordz[n_l,63])/!dtor
 ; ============= Transformacao para o espaco 3D  ======================
 for i=0,n_l-1.+1 do begin
    coord[0,i,*,*]=[[coordx[i,*]],[coordy[i,*]],[coordz[i,*]],[coord1[i,*]]]
;    T3D, /reset,translate=[0.,0.,0.0],Rotate=[0.0,zz[i]+st_sun.a_inc, $
;         90.-st_sun.azim],matrix=tmp
    T3D, /reset, translate=[0.,0.,0.0],Rotate=[0.,zz[i]+st_sun.a_inc, $
         0.],matrix=tmp
    T3D, tmp,Rotate=[st_sun.d_ang,0.0,90.-st_sun.azim],matrix=tmp
    T3D, tmp, translate=[0.0,0.0,+st_sun.L_arc-st_sun.deep],matrix=tmp
    T3D, tmp, Rotate=[-st_sun.b_lat,0.0,0.0],matrix=tmp
    coord[0,i,*,*]=reform(coord[0,i,*,*])#tmp
 endfor
                                ; Coord esta' no sistema de
                                ; coordenadas magnetico. Em 3D no
                                ; meridiano central.
 dip_coord=[[0.],[0.],[0.],[1.0]]      ; Coordenada do dipolo era [[D],
 dip_coord=dip_coord#tmp
; =====================================================================
; for i=0,n_l-1.+1 do begin
;    T3D, /reset, Rotate=[0.0,zz[i],st_sun.azim],matrix=tmp
;    T3D, tmp, translate=[0.0,0.0,st_sun.deep],matrix=tmp
;    coord[0,i,*,*]=reform(coord[0,i,*,*])#tmp
; endfor

; =====================================================================
 ; ====================================================================

 for i=1,n_phase-1 do coord[i,*,*,*]=coord[0,*,*,*] 
                                ; As outras fases nao foram 
                                ; calculadas ainda. Serao logo
                                ; abaixo. Sao uma copia da fase 
                                ; atual giradas.Visao Equatorial.

 ;=================== Rotacao para o eixo de giro do Sol. ===================
 m_phi=st_sun.phi-0.0;st_sun.b_long
 for j=0.,n_phase-1 do begin
    phase=(j*(360.0/n_phase) + m_phi) mod 360.  
    T3D, /reset, Rotate=[0.0,phase,0.],matrix=tmp
    if j eq 0 then dip_coord=dip_coord#tmp
    for i=0,n_l-1+1 do begin
       coord[j,i,*,*]=reform(coord[j,i,*,*])#tmp
    endfor
 endfor
 ;==========================================================
                                ; As coordenadas estao no sistema de
                                ; GIRO e com visao equatorial porem
                                ; com o phi correto.
 ;============================================================================
                                ;       Coordenada do meio da linha
                                ;       central na visao do observador
                                ;       e eixo de giro
 coord_c=reform(coord[0,n_l,*,*])
 T3D, /reset, Rotate=[90.0-st_sun.i,0.0,0.],matrix=tmp 
 coord_c=coord_c#tmp
 x0=coord_c[63,0]   ; pto medio da linha central
 y0=coord_c[63,1] 
 z0=coord_c[63,2] 
 m_r0=sqrt(x0^2+y0^2+z0^2)
 coord=coord[*,0:n_l-1,*,*]
 dip_coord=dip_coord#tmp
 ;===============================Limites do Spot==============================
 xx=fltarr(n_phase,n_l*127)
 yy=fltarr(n_phase,n_l*127)
 zz=fltarr(n_phase,n_l*127)
 x_c=fltarr(n_phase)
 y_c=fltarr(n_phase)
 z_c=fltarr(n_phase)
 x_min=fltarr(n_phase)
 for j=0.,n_phase-1 do begin
    phase=(j*(360.0/n_phase) + m_phi) mod 360.  
    xx[j,*]=reform(coord[j,*,*,0],n_l*127)
    yy[j,*]=reform(coord[j,*,*,1],n_l*127)
    zz[j,*]=reform(coord[j,*,*,2],n_l*127)
    T3D, /reset, Rotate=[0.0,phase-m_phi,0.],matrix=tmp
    T3D, tmp, Rotate=[90.0-st_sun.i,0.0,0.],matrix=tmp 
                                    ; medias feitas na direcao do observador
    tmp_c=[[reform(xx[j,*])],[reform(yy[j,*])],[reform(zz[j,*])],$
                                                     [replicate(1.,n_l*127)]]
;    tmp_c=tmp_c#tmp
    xx[j,*]=tmp_c[*,0]
    yy[j,*]=tmp_c[*,1]
    zz[j,*]=tmp_c[*,2]   
    x_c[j]=mean(xx[j,*]) 
    y_c[j]=mean(yy[j,*])
    z_c[j]=mean(zz[j,*])
;   DEBUG
;    print,'phase=',phase, ' graus'
;    print,'Max.elongacao da fase-> ',max([abs(xx[j,*]-x_c[j]),$
;                abs(yy[j,*]-y_c[j]),abs(zz[j,*]-z_c[j])])
;    mi=min(yy[j,*],max=ma)
;    plot,xx[j,*],yy[j,*],yrange=[mi,ma],xst=1,/iso
;    plots,x_c[j],y_c[j],psy=2
;    plots,x0,y0,psy=2,col=100
;    stop
;   FIM do DEBUG
 endfor  
 for i=0,n_phase-1 do x_min[i]=max([abs(xx[i,*]-x_c[i]),$
                                   abs(yy[i,*]-y_c[i]),abs(zz[i,*]-z_c[i])])
 x_minimo=-max(x_min)*1.1  ;  dez por cento de sobra nas bordas
 center=[[x_c],[y_c],[z_c]]
; ===========================================================================
 ; ===================== Coordenadas so' para a fase atual ==================
 ;    A mascara abaixo nao funciona pois esta' na visao equatorial  
 T3D,/reset,Rotate=[90.0-st_sun.i,0.,0.],matrix=tmp
 for i=0.,n_l-1 do begin 
    cc=reform(coord[0,i,*,*])
    cc=cc#tmp
    dd=reform(sqrt(cc[*,0]^2+cc[*,1]^2))
    tz=reform(cc[*,2])
    pp=where(tz ge 0 or dd ge 1.,compl=npp)
    if npp[0] ne -1 then coord[0,i,npp,*]=!values.f_nan
 endfor
 coord=reform(coord[0,*,*,*])
 cop=coord
; ==============================================================


; ==============================================================
 flag_error=0
 limite={x_min:x_minimo, center:center,coord:cop,x0:x0,y0:y0,z0:z0,$ 
         r0:m_r0,lamb:m_lamb,spot:cc,flag_error:flag_error, $
         dip_coord:dip_coord[0:2]}
end
