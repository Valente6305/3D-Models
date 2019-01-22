pro nobeyama4,dir=dir
; O nome do diretorio dir tem que ter a barra, ex: '~/caca/'

if not keyword_set(dir) then dir='~/HDD2/Modelos_2016/'
restore,'position.sav'
hels=norh.hel
n=n_elements(hels)
ps=strpos(hels,'S')
pn=strpos(hels,'N')
pe=strpos(hels,'E')
pw=strpos(hels,'W')
ips=where(ps ne -1,cs)
ipn=where(pn ne -1,cn)
ipe=where(pe ne -1,ce)
ipw=where(pw ne -1,cw)
flux=reform(norh.flux[0,*])
ss=reverse(sort(flux))
phi=fltarr(n)
b_lat=fltarr(n)

for i=0, cs-1 do b_lat[ips[i]]=-float(strmid(hels[ips[i]],ps[ips[i]]+1,2))
for i=0, cn-1 do b_lat[ipn[i]]=float(strmid(hels[ipn[i]],pn[ipn[i]]+1,2))
for i=0, cw-1 do phi[ipw[i]]=float(strmid(hels[ipw[i]],pw[ipw[i]]+1,2))
for i=0, ce-1 do phi[ipe[i]]=-float(strmid(hels[ipe[i]],pe[ipe[i]]+1,2))

phi=phi[ss]
b_lat=b_lat[ss]

freq=[3.75,9.4,17.0,35.]*1.e9  
fmem=freq
nn=16

;==========Parametros minimos do modelo ===========                                                                                      
 par_min=fltarr(10)
 par_min[0]=0.2 ; ie 
 par_min[1]=1. ; delta       
 par_min[2]=6. ; Ne -> abaixo de 10^6 nao deve ser importante                                                                                  
 par_min[3]=-50. ; Assimetria em percentis                                                                   
 par_min[4]=-90. ; Azimuth Graus                                                                                                               
 par_min[5]=.003 ; Raio do Apex                                                                                                               
 par_min[6]=0.008 ; Foot separation                                                                                                            
 par_min[7]=0.008 ; h_arc                                                                                                                      
 par_min[8]=-40. ; Inclinacao                                                                                                                  
 par_min[9]=800. ; Bmag                                                                                                                        
                            
;======= Parametros maximos do modelo ===========                                                                                    
 par_max=fltarr(10)
 par_max[0]=2. ; ie                                                                                                                            
 par_max[1]=4.4 ; delta                                                                                                                                                
 par_max[2]=8.  ; Ne -> abaixo de 10^6 nao deve ser importante                                                                                  
 par_max[3]=50.                 ; Assimetria em percentis    
 par_max[4]=90. ; Azimuth Graus                                                                                                                
 par_max[5]=.02 ; Raio do Apex                                                                                                                
 par_max[6]=0.06                ; Foot separation                                                                                                             
 par_max[7]=0.06                ; h_arc                                                                                                                       
 par_max[8]=40.                 ; Inclinacao                                                                                                                  
 par_max[9]=3500.               ; Bmag                                                                                                                       
                           
 
energy=[10.,1.e5]
code=1
mm=0
T=1.e8
ip=[0,1,2,3,4]
;print,'frequencia=', freq/1.e9

temporario=1
while temporario do begin
   denovo:
;print,'denovo'
   i=randomu(seed)*(n-1l)
   in0=randomu(seed)*(par_max[8]-par_min[8])+par_min[8]
   as0=(randomu(seed)*(par_max[3]-par_min[3])+par_min[3])/100.
   mag=randomu(seed)*(par_max[9]-par_min[9])+par_min[9]
   delta=randomu(seed)*(par_max[1]-par_min[1])+par_min[1]
   nel=10.^(randomu(seed)*(par_max[2]-par_min[2])+par_min[2])
   phi0=phi[i]
   b_lat0=b_lat[i]
   np0=fix(randomu(seed)*4.999)
   azim0=randomu(seed)*(par_max[4]-par_min[4])+par_min[4]
   r_arc0=randomu(seed)*(par_max[5]-par_min[5])+par_min[5]
   foot_s0=randomu(seed)*(par_max[6]-par_min[6])+par_min[6]
   h_arc0=randomu(seed)*(par_max[7]-par_min[7])+par_min[7]
   ih=randomu(seed)*(par_max[0]-par_min[0])+par_min[0]
   npi=ip[np0]
   dens=10.^(randomu(seed)*(12.-9.)+9.)
   ie0=ih
   mostra=0
   
;   if mag le 2700. then begin
;      goto, denovo   
  ; endif

  
   if r_arc0 ge h_arc0 then begin
      print,'limite r_arc'
      goto, denovo   
   endif
   phi_limit=atan(r_arc0/h_arc0)/!dtor
   if ((phi_limit+abs(in0))*1.05) gt 90. then begin
      print,'limite_phi'
      goto,denovo 
   endif
        
   seed_file,st, /new, nn=nn, phi=phi0,b_lat=b_lat0,in=in0,assim=as0, $
              ie=ie0,ip=npi,r_arc=r_arc0,foot_s=foot_s0,np=dens $
              ,h_arc=h_arc0,azim=azim0,nel=nel,t=t,energy=energy, $
              delta=delta, mm=mm, mag=mag, code=code, freq=freq

   ofile=loop_file(st=st)
   freq=fmem
   if file_search(dir+ofile) eq '' then begin
            calc_image,imagem,st=st, f_save=f_save, mostra=mostra, np=np $
               ,px2cm=px2cm, mag_res=mag_res, ofile=ofile, ind=npi, dir=dir
            print,ofile 
         endif else print,'Ja existe',' ',ofile
pula:
     
endwhile

end
 
