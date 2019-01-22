function fill_np,mag_res,high
time0=systime(/sec)
  if total(size(high)) eq 0 then high=0
  n=mag_res.n
  st_arc=mag_res.st_arc
  lim=mag_res.lim
  r_v=fltarr(n,n,n)
  n_p=dblarr(n,n,n)
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
  T3D,/Reset,Translate=reform(ccc[0,0:2]),matrix=tmp0
  i_z=0l
  for t_z=spc_min,spc_max, delta do begin
    c=[[reform(t_x,n^2)],[reform(t_y,n^2)],[replicate(t_z,n^2)],$
     [replicate(1.0,n^2)]]   ; Coordenadas em pixels
    c=c#tmp0                 ; Aqui coord.s no lugar espacial do momento
                             ; com origem no centro da estrela (Sol)
    rvc=sqrt(c[*,0]^2+c[*,1]^2+c[*,2]^2)/st_arc.r_sun ;raio vetor c/origem centro
    r_v[*,*,i_z]=reform(rvc,n,n)  ; Raio vetor do pto (origem no centro)
    i_z=i_z+1l
 endfor
 pp=where(reform(mag_res.b) ne 0)
 if pp[0] ne -1 then r_v[pp]=r_v[pp] else return,-1 
 ;Precisa fazer um case abaixo
SWITCH high OF 
   0: BEGIN 
          n0=1.d10 
          h0=0.00434294    
          BREAK 
      END 
   1: BEGIN 
          n0=1.d11
          h0=0.00206418      
          BREAK 
      END 
   2: BEGIN 
          n0=1.d12 
          h0=0.00289530      
          BREAK 
      END 
   3: BEGIN 
          n0=10.*1.d12 
          h0=0.00289530      
          BREAK 
      END
   4: BEGIN 
          n0=20.d12 
          h0=0.00289530      
          BREAK 
      END 
ELSE: BEGIN 
          n0=1.d10 
          h0=0.00434294    
          BREAK 
      END 
ENDSWITCH 


; if high eq 3 then begin
;    n0=2.*5.*1.d12 ; o cinco e' para o ind=2, o 2 e' ind=3
;    h0=0.00289530
; endif else begin
;    n0=1.d10
;    h0=0.00434294
; endelse 
n_p[pp]=n0*exp(-(r_v[pp]-st_arc.l_arc)/h0)
timef=systime(/sec)-time0    
return,n_p

end
