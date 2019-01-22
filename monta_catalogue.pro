pro monta_catalogue

f=file_search('~/HDD2/Modelos_2016/*.sav'); Diretorio dos modelos       
ff=f[0,50000]
file_distrib, dir, st, names=ff
pp03=where(st.ie le 0.5)
pp12=where(st.ie gt 0.5 and st.ie le 1.5)
pp19=where(st.ie gt 1.5 )

st03=st[pp03]
st12=st[pp12]
st19=st[pp19]

ff03=ff[pp03]
ff12=ff[pp12]
ff19=ff[pp19]

ss03=sort(st03.phi)
ss12=sort(st12.phi)
ss19=sort(st19.phi)

ff03=ff03[ss03]
ff12=ff12[ss12]
ff19=ff19[ss19]

freq=[3.75,9.4,17.0,35.]
names=[ff03,ff12,ff19]
n_mod=[n_elements(ff03),n_elements(ff12),n_elements(ff19)]
file_distrib, dir, st, names=names
print,'Numero de modelos=',n_elements(st)
flux_distr,st,names,lat,long,flux,area,w,file=file ; Demora 
n=n_elements(lat)
t=st[0]
tt={n:t.n,phi:t.phi,B_lat:t.b_lat,Az:t.az,in:t.in,as:t.as,ie:t.ie,ip:t.ip,r_arc:t.r_arc,$
          foot_s:t.foot_s,h_arc:t.h_arc,Nel:t.Nel,T:t.t,En:t.en, $
          delta:t.delta,m:t.m,B:t.b,code:t.code, freq:'',image:fltarr(32,32,4),$
          flux:fltarr(4),area:fltarr(6),file:''}
st_models=replicate(tt,n)
for j=0l,n-1 do begin 
   for i=0,17 do begin 
      st_models[j].(i)=st[j].(i)
   endfor

   st_models[j].image=reform(w[*,*,*,j])
   st_models[j].flux=reform(flux[*,j])
   st_models[j].area=reform(area[*,j])
   st_models[j].file=names[j]
endfor
save,file='Kawate_catalogue_2016.sav',n_mod,st_models
end


