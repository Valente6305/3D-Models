Function jsgn,x
  ; Funcao que pega o sinal do vetor x
  pp=where(x eq 0.,compl=np)
  sgn=x*0.0+1.0
  if np[0] ne -1 then sgn[np]=abs(x[np])/float(x[np])
return,sgn
end
;function dis,x,par=par, max=max
;    if keyword_set(max) then fact=max else fact=1.0
;   return,1.5*(10.125*x^4+par)*exp(-4.5*x^2)/(0.41873+1.2449*par)
;   return,(2.0*(x/fact)^4+par)*exp(-2.0*(x/fact)^2)/(0.41873+1.2449*par)
;end

FUNCTION dis,x,par=par,max=max,bal=bal
        if not keyword_set(bal) then ampl=1.0 else ampl=bal
        if keyword_set(max) then fact=max else fact=1.0
        if par gt 1 then begin
           x0=(1.*par-1.0)<1.0
           den=(2.0-par)>0.02
        endif else begin
           x0=0.0
           den=par>0.02
        endelse
        tmp=ampl*exp(-(x-1.0+x0)^2/den)+exp(-(x+1.0-x0)^2/den)
        tmp=tmp/max(tmp)*fact
        return,tmp
END

pro mag_sun,st_arc,n,b,theta,comp,px,mag_res=mag_res, giro=giro,nel=nel,T=T, $
              var_t=var_t,var_n=var_n,error=error,lim=lim,r_v=r_v, n_p=dens_p, $
              debug=debug, par=par, leng=leng


; J E R Costa - 01Jan2012
; Rotina que calcula os campos magneticos nas coordenadas de um cubo 
; de voxels centrado no centro ode giro do loop. 
; Os voxels no cubo tem px diferente definido no st_sun 
; (default) do raio da Sol. O px que sai aqui corrige isso.
  time0=systime(/sec)      
  if not keyword_set(nel) then nel=2.e7
  if not keyword_set(T) then T=1.e8 ; Kelvin (ambiente dentro do loop)
  if not keyword_set(var_t) then var_t=0 ; Temp. variavel no espaco (N.Impl.Y.)
  if not keyword_set(var_n) then var_n=0 ; Dens. amb. variavel no espaco (N.Impl.Y.)

; tornando n uma variável long para evitar overflow de n*n
  error=0
  n=long(n)
  dens=nel
  temper=T
  ;frac_r=st_arc.frac_r
; Determinando os limites do arco
  tmp_st=st_arc
;-----------------------------------------------------------------------
;============== Determinacao dos LIMITES do ARCO =======================
  limite_arc,tmp_st,lim,n_linhas=n
  if lim.flag_error then begin 
     error=1
     goto, sem_limite
  endif
;-----------------------------------------------------------------------
  tem=0
; Criando a variavel coord dos voxels em torno do centro 
; do Obj (tubo)
  leng=fltarr(n,n,n)
  tleng=fltarr(n*n)
  r_v=fltarr(n,n,n)
  b=fltarr(n,n,n)
  nel=fltarr(n,n,n)
  T=fltarr(n,n,n)
  comp=fltarr(n,n,n)
  giro=fltarr(n,n,n)
  theta=fltarr(n,n,n)
  n_p=fltarr(n,n,n)
  delta=2.*abs(lim.x_min)/(n-1.)*st_arc.r_sun ; Esta' em pixel
  px=delta*st_arc.px
  spc_min=-(n-1-(n-1)/2.)*delta
  spc_max=(n-1-(n-1)/2.)*delta
  t_x=((findgen(n)-(n-1.)/2.)#replicate(1,n))*delta
  t_y=replicate(1,n)#(findgen(n)-(n-1)/2.)*delta

  ; Obtencao do centro
  ccc=reform(lim.center[0,0:2]*st_arc.r_sun)   ; centro do cubo do loop. Visao do obs
  ccc=[[ccc[0]],[ccc[1]],[ccc[2]],[1.]];       ; Visao do observador (checado) 
                                               ; no eixo de giro (em r_sun)
  ; Testei o pto central da meia altura do arco = ficou pior!
;  rrrr=1.+st_arc.h_arc/2.
;  xxxx=rrrr*cos(!dtor*st_arc.b_lat)*sin(!dtor*st_arc.phi)*st_arc.r_sun
;  yyyy=rrrr*sin(!dtor*st_arc.b_lat)*st_arc.r_sun
;  zzzz=rrrr*cos(!dtor*st_arc.b_lat)*cos(!dtor*st_arc.phi)*st_arc.r_sun
;  ccc=[[xxxx],[yyyy],[zzzz],[1.]]
;stop
  T3D,/Reset,Translate=reform(ccc[0,0:2]),matrix=tmp0

;  Qual voxel representa o centro de referencia?
  dx=round((lim.x0*st_arc.r_sun-ccc[0,0])/delta+(n-1.)/2.)
  dy=round((lim.y0*st_arc.r_sun-ccc[0,1])/delta+(n-1.)/2.)
  dz=round((lim.z0*st_arc.r_sun-ccc[0,2])/delta+(n-1.)/2.)
  ref={i:dx,j:dy,k:dz}
; Constante Magnetica 
;  B0=(lim.r0*st_arc.r_sun)^3/sqrt(1.0+3.0*(sin(lim.lamb*!dtor))^2)*st_arc.B
;  dipole=(st_arc.deep+st_arc.L_arc-1.0)*st_arc.r_sun
;
  dipole=sqrt((st_arc.deep-st_arc.l_arc)^2+st_arc.l_arc^2-$
         2*sqrt(1.0-(st_arc.foot_s/(2.*st_arc.l_arc))^2)* $
              st_arc.l_arc*(st_arc.l_arc-st_arc.deep))*st_arc.r_sun
 B0=(dipole)^3/2.d0*st_arc.B
 ;B0=st_arc.b0
 i_z=0l

  ;B0=1.0
  ; Comprimento da linha central
   elonga=(st_arc.deep+st_arc.h_arc)
   ang=atan(st_arc.foot_s/2.,st_arc.deep)*findgen(20)/19.0
   rrr=elonga*(cos(ang))^2
   xxx=rrr*cos(ang)
   yyy=rrr*sin(ang)
   compmax=((compr(xxx,yyy))[19])*st_arc.r_sun
;
;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@INICIO DO LOOP@@@@@@@@@@@@@@@@@@@@@@@@@
  for t_z=spc_min,spc_max, delta do begin

  c=[[reform(t_x,n^2)],[reform(t_y,n^2)],[replicate(t_z,n^2)],$
     [replicate(1.0,n^2)]]   ; Coordenadas em pixels
  c=c#tmp0                   ; Aqui coord.s no lugar espacial do momento
                             ; com origem no centro da estrela (Sol)
  rvc=sqrt(c[*,0]^2+c[*,1]^2+c[*,2]^2)/st_arc.r_sun ;raio vetor c/origem centro
 
  eclipse=((sqrt(c[*,0]^2+c[*,1]^2) ge st_arc.r_sun) or (c[*,2] ge 0.))
;-----------------------------------------------------------------------
; As coordenadas do cubo (c) ate' aqui estao corretas para o primeiro Z
; e o centro fica realmente no centro do cubo. Se o centro estiver 
; correto o resto tbem esta'. O cubo tem seu eixo Z apontado para 
; o observador.Obvio o cubo e' um objeto 3D.
;-----------------------------------------------------------------------

; Centro do disco plano no limbo Oeste (x0,y0,z0)
; O Anchor esta' no espaco e no limbo oeste. Portanto, o anchor[0,*]
; esta' bem no limbo magnetico.
  rho=st_arc.l_arc*st_arc.r_sun          ; Altura da base do arco em pixel
  rs=st_arc.r_arc*st_arc.r_sun           ; raio do arco em pixel
  rss=(st_arc.l_arc+st_arc.h_arc)*st_arc.r_sun ; Altura do topo do arco
;---------------------------------------------------------------------------

; Giro solido para visao equatorial (visao equator. mesmo)
;  T3D,/Reset,Rotate=[arc.st.i-90.,0.,0.],matrix=tmp 
; -(arc.st.phi) zera o phi (meridiano central)
  T3D,/Reset,Rotate=[st_arc.i-90.,-st_arc.phi,0.],matrix=tmp 
  c=c#tmp ; visao equatorial e no meridiano Z. Ou seja o centro do plano 
;           vem para o meridiano do meio dia. Fica no entorno.
  c=c#invert(st_arc.tr) ; Esta' no eixo magnetico
;Estas coords c desta fatia tem origem no dipolo e o dipolo 
; esta' vertical(y) e com centro no equador 
;---------------------------------------------------------------------------
; A operacao acima coloca o cubo (fatia de 3D para o Z do momento) com 
; visao equatorial e no eixo magnetico. Porem, nao voltou para phi!?!
;---------------------------------------------------------------------------
;---------------------------------------------------------------------------
; CALCULO DA EQUACAO DAS LINHAS MAGNETICAS
  r=sqrt(c[*,0]^2+c[*,1]^2+c[*,2]^2)       ; raio vetor de cada voxel
;  rd=sqrt((c[*,0]-lim.dip_coord[0])^2+(c[*,1]-lim.dip_coord[1])^2+(c[*,2]-lim.dip_coord[2])^2)
  lamb=atan(c[*,1]/sqrt(c[*,0]^2+c[*,2]^2)); lat mag. de cada voxel  
  lon=!dpi/2.-atan(c[*,0],c[*,2])          ; qto falta pro limbo Oeste.
  r0=r/((cos(lamb))^2)                     ; ro da linha de cada pto
  tb=b0*sqrt(1.0+3.0*(sin(lamb))^2)/r^3    ; Campo Magnetico
  tr_v=r

;  lamb0=acos(sqrt((st_arc.l_arc-1.0+st_arc.deep)*st_arc.r_sun/r0)) ; lamb do pe' do arco
  lamb0=acos(sqrt((st_arc.deep)*st_arc.r_sun/r0)) ; lamb do pe' do arco
 ; if i_z eq 0 then begin 
 ;    b0=st_arc.b/max(tb)
 ;    tb=b0*sqrt(1.0+3.0*(sin(lamb))^2)/r^3
 ;    print,'B0=',B0
 ; endif
;  tb=r*0.0+1.0
  tdens=tb*0.0+dens
  ttemp=tb*0.0+temper
  tn_p=tb*0.0+dens_p
 

  ; Determinando ptos internos do arco
  dss=sqrt((c[*,0]-0)^2+(c[*,1]-0)^2+(c[*,2]-0)^2) ; raio vetor do pto(pixel)

  raio=st_arc.r_arc*st_arc.r_sun  ; em pixel
  z=r0*cos(!dpi/2.-lon)  ; em pixel
  y=z*0.0  ; em pixel
  x=r0*sin(!dpi/2.-lon)  ; em pixel
; ==== lim.xyz0  sao as coords do APEX da lonha central ============
  cen=[[lim.x0],[lim.y0],[lim.z0],[1.0]] ; Esta' em relacao ao centro do 
                                ; sol e precisa estar em relacao ao
                                ; centro do dipolo
  cen[*,0:2]=cen[*,0:2]*st_arc.r_sun

  cen=cen#tmp ; visao equatorial e no meridiano Z
  cen=cen#invert(st_arc.tr) ; Esta' no eixo magnetico
  x0=cen[0,0];*st_arc.r_sun
  y0=cen[0,1];*st_arc.r_sun
  z0=cen[0,2];*st_arc.r_sun

  ds=sqrt((x-x0)^2+(y-y0)^2+(z-z0)^2)  ; em pixel
;  ==========  Grafico para debug =================================
  if keyword_set(debug) then begin
     if i_z eq 0 then v_ds=fltarr(n,n,n)
     if i_z eq 0 then v_la=fltarr(n,n,n)
     if i_z eq 0 then window,0,xsize=600,ysize=500
     era0=bytarr(6*25,8)
     ww=((ds le raio) and (eclipse) and (rvc ge st_arc.l_arc));(tr_v ge dipole))
     if i_z eq 0 then plot,t_x,t_y,/iso,psy=3
     pp_w=where(ww eq 1)
     t_tb=reform(tb,n,n)
     if pp_w[0] ne -1 then begin
        t_pp=where(t_tb[pp_w] gt st_arc.b*(1.+st_arc.assim))
        plots,t_x[pp_w],t_y[pp_w],psy=7,col=100 
        plots,t_x[pp_w],t_y[pp_w],psy=6,col=100 
        if t_pp[0] ne -1 then begin
           print,string(7B),'Campo alto'
           plots,t_x[pp_w[t_pp]],t_y[pp_w[t_pp]],psy=7,col=150 
           plots,t_x[pp_w[t_pp]],t_y[pp_w[t_pp]],psy=6,col=150 
        endif
        v_ds[*,*,i_z]=reform(r,n,n)  
        v_la[*,*,i_z]=reform(lamb,n,n)  
        tv,era0,0.5*!d.x_size,0.97*!d.y_size
        xyouts,.5,0.97,'ds_max='+string(max(v_ds[*,*,i_z])),/norm
        tv,era0,0.1*!d.x_size,0.97*!d.y_size
        xyouts,.1,0.97,'ds_min='+string(min(v_ds[*,*,i_z])),/norm
        oo=where(v_ds ne 0)
        if oo[0] ne -1 then dsmin=min(v_ds[oo],max=mads,imin,subscript_max=imax)
        tv,era0,0.5*!d.x_size,0.02*!d.y_size
        xyouts,.5,0.02,'ds_max='+string(mads),/norm
        tv,era0,0.1*!d.x_size,0.02*!d.y_size
        xyouts,.1,0.02,'ds_min='+string(dsmin),/norm
        tv,era0,0.1*!d.x_size,0.08*!d.y_size
        xyouts,.1,0.08,'B_min='+string(b0*sqrt(1.0+3.0*(sin(v_la[oo[imin]]))^2)/v_ds[oo[imin]]^3 ),/norm
        wait,0.3
     endif
  endif
;    if i_z eq 16 then stop
;===================================================================
;  REGION OF INTEREST
;  pp=where((ds le raio) and (eclipse) and (tr_v ge dipole) , compl=np)
  pp=where((ds le raio) and (eclipse) and (rvc ge st_arc.l_arc) , compl=np)
 
; (ds le raio) -> linhas com distancia de r0 ate r0 central < r_arc
; (dss le rss) -> pontos abaixo dos pes do loop (nao usado!)
; (eclipse)    -> pontos fora do eclipse pelo Sol
; (tr_v ge dipole) -> ptos acima do raio vetor do pé do arco
; (rvc ge st_arc.l_arc)-> ptos acima do raio vetor do pé do arco
;===================================================================
  if np[0] ne -1 then begin 
     tb[np]=0.0 
     tdens[np]=0.0 
     ttemp[np]=0.0 
     tr_v[np]=0.0
     tn_p[np]=0.0
  endif
  b[*,*,i_z]=reform(tb,n,n)           ; Campos magneticos
  nel[*,*,i_z]=reform(tdens,n,n)      ; Densidade
  T[*,*,i_z]=reform(ttemp,n,n)        ; Temperatura
;  r_v[*,*,i_z]=reform(tr_v,n,n)      ; Raio vetor do pto (origem no dipolo)
  r_v[*,*,i_z]=reform(rvc,n,n)        ; Raio vetor do pto (origem no centro)
  n_p[*,*,i_z]=reform(tn_p,n,n)       ; Densidade ambiente de protons

;================== Fim do Calculo de B ============================

;Componentes dos versores
vx=double(tb*0.0)
vy=vx
vz=vx

;=============== Inicio do Calculo de theta ========================
; Direcao do Campo (tangente da linha no Limbo Oeste)
; Girando todos (longitudes) para o plano de referencia (plano Oeste)
;
  us=20 ; Numero de coord. na linha para calcular o tamanho dela
  ang_f=findgen(us+1)/float(us) ; positivo pq vai so' ate' apex
  compt=fltarr(n_elements(r0))
;  comp_l=compt
  for i=0,n_elements(r0)-1l do begin ; Leva C para o limbo oeste
     ; Comprimento maximo de cada linha comp_l
     ang=lamb0[i]*findgen(20)/19.0
     rrr=r0[i]*(cos(ang))^2
     xxx=rrr*cos(ang)
     yyy=rrr*sin(ang)
     tleng[i]=((compr(xxx,yyy))[19]) ; *st_arc.r_sun

     T3D,/Reset,Rotate=[0.,lon[i]/!dtor,0.],matrix=tmp
     c[i,*]=reform(c[i,*])#tmp
     xx_l=r0[i]*(cos(ang_f*abs(lamb[i])))^3
     yy_l=r0[i]*(cos(ang_f*abs(lamb[i])))^2*sin(ang_f*abs(lamb[i]))
     compt[i]=(compr(xx_l,yy_l,px=px,mm=tleng[i]))[us]
 ; if finite(compt[i]) eq 0 then stop
  endfor 

  if np[0] ne -1 then compt[np]=0.
  comp[*,*,i_z]=reform(compt,n,n)       ; comprimentos ate o apex
  if np[0] ne -1 then tleng[np]=0.0
  leng[*,*,i_z]=reform(tleng,n,n)

  ; Testei: lamb e' exatamente igual atan(c[*,1],c[*,0])
  rr0=double(r0)
  theta0=double(lamb)   ; Angulo polar de cada pto no Limbo
  tp=where(abs(theta0) lt 1e-4,compl=ntp) ; linha nova !!!!!!
  ta=double(theta0*0.0)
  if tp[0] ne -1 then $
      ta[tp]=-1.d0/theta0[tp]  ; Da'um numero gde (B aponta p/baixo -/+)
  if ntp[0] ne -1 then begin
        ;  Linha tangente
        q=(3.d0*c[ntp,0]^(1./3)*((c[ntp,0])^(2./3)-rr0[ntp]^(2./3)))
        pq=where(theta0[ntp] lt 0.)
        if pq[0] ne -1 then q[pq]=-q[pq]
        pq=where(q ne 0.,compl=npq)
        if pq[0] ne -1 then ta[ntp[pq]]=rr0[ntp[pq]]^(1./3)* $
                    (3.*(c[ntp[pq],0])^(2./3)-2.*rr0[ntp[pq]]^(2./3))* $
                    sqrt((rr0[ntp[pq]]^(2./3)-(c[ntp[pq],0])^(2./3))/$
                    (rr0[ntp[pq]]^(2./3)))/q[pq] 
        if npq[0] ne -1 then ta[ntp[npq]]=1.e7
  endif  
ttheta=theta0*0.0
tgiro=ttheta          ; Dimensao do tgiro com zero onde x e y = 0

if pp[0] ne -1 then begin
   vx[pp]=cos(atan(ta[pp]))*jsgn(reform(c[pp,1]))
   vy[pp]=sin(atan(ta[pp]))*jsgn(reform(c[pp,1]))
   vc=[[0.,vx[pp]],[0.,vy[pp]],[0.,vz[pp]],[replicate(1.0,n_elements(pp)+1)]]
vcc=vc[1:*,*]
tmpt=fltarr(n_elements(pp),4,4)
tpp=[0,pp]
tlon=[0.0,lon]
;   for i=0l,n_elements(pp)-1l do begin 
;      T3D,/Reset,Rotate=[0.,-lon[pp[i]]/!dtor,0.],matrix=tmp
;      T3D,tmp#st_arc.tr,Rotate=[0.,st_arc.phi,0.],matrix=tmp
;      vc[i,*]=reform(double(vc[i,*]))#tmp
;      if i ne 0 then vc[i,0:2]=vc[i,0:2]-vc[0,0:2]
;   endfor 
   for i=0l,n_elements(pp) do begin 
      T3D,/Reset,Rotate=[0.,-tlon[tpp[i]]/!dtor,0.],matrix=tmp
      T3D,tmp#st_arc.tr,Rotate=[0.,st_arc.phi,0.],matrix=tmp
      vc[i,*]=reform(double(vc[i,*]))#tmp
      if i ne 0 then vc[i,0:2]=vc[i,0:2]-vc[0,0:2]
   endfor 
   vc=vc[1:*,*]
   T3D,/Reset,Translate=[0.0,0.0,0.0],matrix=tmp
   T3D,tmp,Rotate=[90.0-st_arc.i,0.,0.],matrix=tmp ; sempre direcao do obs.
   vc=vc#tmp
   ; Fazendo o Calculo do angulo de visada
   cs=acos(vc[*,2]) ; pelo produto escalar
   ttheta[pp]=cs
   tgiro[pp]=(atan(-vc[*,0],vc[*,1])+2*!dpi) mod (2.*!dpi);In radianos
endif 

theta[*,*,i_z]=reform(ttheta,n,n)
giro[*,*,i_z]=reform(tgiro,n,n) ; Para calcular o angulo de polarizacao
i_z=i_z+1l
;================= Fim do Calculo de theta ========================
endfor 

;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ 
b_nz=where(b ne 0)
if b_nz[0] eq -1 then error=1 else begin 
   ; Correcao de B para o maximo pedido = st_arc.b
   bmax=max(b)
   b[b_nz]=b[b_nz]*st_arc.b*(1.+abs(st_arc.assim))/bmax
   print,'Novo maximo de B=',max(b)
   print,'B0 mudou na estrutura do arco!!!  Era:', st_arc.b0
   st_arc.b0=st_arc.b0*st_arc.b*(1.+abs(st_arc.assim))/bmax
   print, 'Ficou:',st_arc.b0
   ; Para garantir que ref existe dentro do volume de voxels
   if b[ref.i,ref.j,ref.k] eq 0 then begin
      print, 'Ref mudou de i:',ref.i,' j:',ref.j,' k:',ref.k, ' para'
      ind=float(array_indices(b,b_nz))
      d=sqrt((ind[0,*]-ref.i)^2+(ind[1,*]-ref.j)^2+(ind[2,*]-ref.k)^2)
      tmp=min(d,ii)
      ref.i=ind[0,ii]
      ref.j=ind[1,ii]
      ref.k=ind[2,ii]
      print, '             i:',ref.i,' j:',ref.j,' k:',ref.k 
   endif
endelse
sem_limite:
if not error then begin
   peso=dis(comp[b_nz],par=par)
   nel[b_nz]=dens*peso
   fact=n_elements(b_nz)*dens/total(nel[b_nz])
   nel[b_nz]=fact*nel[b_nz]
endif
mag_res={n:n,st_arc:st_arc,b:b,theta:theta/!dtor, comp:comp,giro:giro,T:T,nel:nel,np:n_p, $
         px:px, var_t:var_t, var_n:var_n,lim:lim,raio_v:r_v,line_len:compmax}

timef=systime(/sec)-time0
;print,'Tempo do mag_sun: ', timef 

end


