pro converg_Y

;======= Parametros minimos do modelo ===========
;================================================                                    
 par_min=fltarr(10)
 par_min[0]=0.2 ; ie                                                   
 par_min[1]=1. ; delta                                                               
 par_min[2]=10.^6. ; Ne -> abaixo de 10^6 nao deve ser importante                        
 ;par_min[3]=7. ; T  -> 10^6 K um loop frio                                           
 ;par_min[4]=6. ; Np -> acima de 10^8 cm-3                       
 par_min[3]=-50 ; Assimetria em percentis
 par_min[4]=-45. ; Azimuth Graus      
 par_min[5]=.002  ; Raio do Apex
 par_min[6]=0.008 ; Foot separation
 par_min[7]=0.008 ; h_arc                                                            
 par_min[8]=-40. ; Inclinacao
 par_min[9]=800. ; Bmag  
 ;par_min[12]=-90. ; Phi 
 ;par_min[13]=-90. ; B_lat        

;======= Parametros maximos do modelo ===========
;================================================                                    
 par_max=fltarr(10)
 par_max[0]=2. ; ie                                                   
 par_max[1]=4.4 ; delta                                                             
 par_max[2]=10.^8. ; Ne -> abaixo de 10^6 nao deve ser importante                        
 ;par_max[3]=9. ; T  -> 10^6 K um loop frio                                           
 ;par_max[4]=12. ; Np -> acima de 10^8 cm-3                       
 par_max[3]=50 ; Assimetria em percentis
 par_max[4]=45. ; Azimuth Graus      
 par_max[5]=.02  ; Raio do Apex
 par_max[6]=0.06 ; Foot separation
 par_max[7]=0.06 ; h_arc                                                             
 par_max[8]=40. ; Inclinacao
 par_max[9]=3500. ; Bmag  
 ;par_max[12]=90.   ; Phi 
 ;par_max[13]=90.  ; B_lat        

;close,1
;openw,1,'YY_parametros_44668.txt'
;printf,1,'ie,delta,Nel,Assim,Azim,Apex_radii ,Foot_s,Height,Incl, Magnetic' 
 restore,'Resultado_136_banco_163mil_16x16.sav'
n=n_elements(struttotal.ofile)
dir1='~/Modelos_teste/'

par_a=fltarr(10,n)
par_o=fltarr(10,n)
VYY=fltarr(10,n)
YY=fltarr(n)

for i=0,n-1 do begin
 
   restore,dir1+struttotal[i].ofile
   pp=strpos(ofile,'nn')
   Nel=10.^(float(strmid(ofile,pp+53,3))/100.)


   if (where(tag_names(mag_res.st_arc) eq 'ASSIM'))[0] eq -1 then tass=0.0 else $
              tass=100.*mag_res.st_arc.assim
   
   par_a[*,i]=[par,delta,Nel,tass,mag_res.st_arc.azim, mag_res.st_arc.r_arc,mag_res.st_arc.foot_s,mag_res.st_arc.h_arc, mag_res.st_arc.a_inc,mag_res.st_arc.b ]
   
   
                                ;for j=0,0 do begin
   restore,struttotal[i].file[0]
   pp=strpos(ofile,'nn')
   Nel=10.^(float(strmid(ofile,pp+53,3))/100.)
   if (where(tag_names(mag_res.st_arc) eq 'ASSIM'))[0] eq -1 then tass=0.0 else $
      tass=100.*mag_res.st_arc.assim
   
   par_o[*,i]=[par,delta,mag_res.st_arc.Nel,tass,mag_res.st_arc.azim,mag_res.st_arc.r_arc,mag_res.st_arc.foot_s,mag_res.st_arc.h_arc, mag_res.st_arc.a_inc,mag_res.st_arc.b ]
   
   VYY[*,i]=(par_a[*,i]-par_o[*,i])/(par_max-par_min) ; guardar
   
   YY[i]=total(abs(VYY[*,i]))

   ;print,'par_a:',par_a
   ;print,'par_o:',par_o
   ;print,'vyy: ', abs(VYY)
   ;print,'YY: ',i, yy
   ;printf,1,i,par_a[*,i]
   ;printf,1,i,par_o[*,i]
   ;printf,1,i,abs(VYY[*,i])
                                ;print,VYY ,YY      
                                ; endfor
endfor
close,1
save,filename='melhor_y_136_chi_163mil_16x16.sav',VYY,YY,par_a,par_o

end
