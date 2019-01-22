pro fits

restore,'evento_2013FEB17.save'
;Calculo das inclinaçoes dos espectros
a1=alog(PICO_FLUX_8_8GHZ[1]/PICO_FLUX_5GHZ[0]); fluxos
a2=alog(fq[1]/fq[0]);frq
alfa=a1/a2 ;radio spectral index
y=[PICO_FLUX_8_8GHZ[1],PICO_FLUX_5GHZ[0]]
K=y[0]/fq[1]^alfa; constante de proporcionalidade
F3=k*3.7^alfa; Fluxo extrapolado
print,'Alfa em 3.7GHz:',alfa ,'   ','K:',K,'   ','Fluxo extrapolado em 3.7GHz:',F3

a1=alog(PICO_FLUX_15GHZ[2]/PICO_FLUX_8_8GHZ[1]); fluxos
a2=alog(fq[2]/fq[1]);frq
alfa=a1/a2 ;radio spectral index
y=[PICO_FLUX_15GHZ[2],PICO_FLUX_8_8GHZ[1]]
K=y[0]/fq[2]^alfa; constante de normalização
F9=k*9.4^alfa; Fluxo extrapolado
print,'Alfa em 9.4GHz:',alfa ,'   ','K:',K,'   ','Fluxo extrapolado em 9.4GHz:',F9

a1=alog(PICO_FLUX_45GHZ[3]/PICO_FLUX_15GHZ[2]); fluxos
a2=alog(fq[3]/fq[2]);frq
alfa=a1/a2 ;radio spectral index
y=[PICO_FLUX_45GHZ[3],PICO_FLUX_15GHZ[2]]
K=y[0]/fq[3]^alfa; constante de normalização
F17=k*17.0^alfa; Fluxo extrapolado
print,'Alfa em 17GHz:',alfa ,'   ','K:',K,'   ','Fluxo extrapolado em 17GHz:',F17

F34=k*34.0^alfa; Fluxo extrapolado
print,'Alfa em 34GHz:',alfa ,'   ','K:',K,'   ','Fluxo extrapolado em 34GHz:',F34

s='N12E20'

tt={freq:fltarr(4),oflux:fltarr(4),flux:fltarr(4,12),sigma:fltarr(12),image:fltarr(32,32,4,12),posic:'',file:strarr(12)}
struttotal=replicate(tt,1)


   search_catalogue,[F3,F9,F17,F34],s,res,catal=catal,ii=ind
   for j=0,11 do struttotal.image[*,*,*,j]=catal[ind[j]].image
   struttotal.freq=res.freq
   for j=0,11 do struttotal.file[j]=catal[ind[j]].file
   struttotal.oflux=reform([F3,F9,F17,F34])
   struttotal.flux=res.flux
   struttotal.sigma=res.sigma
   struttotal.posic=s

save,file='Catalogo_2013FEB17.sav',struttotal

end



