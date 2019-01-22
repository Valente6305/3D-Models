pro roda_calc

freq=[3.75,9.4,17.0,35.0]*1.e9
fmem=freq
h=[0.05,0.34,2.0]
nf=n_elements(freq)
fmt='('+string(nf)+'i4)'
w=string(freq/1e8,format=fmt)
bw=byte(w)
pp=where(bw eq 32)
if pp[0] ne -1 then bw[pp]=48
w='fr'+string(bw)
phiv=[10.0,20.0,-10.0,-20.0]
azimv=[30.0,10.0,-10.0]
b_latv=[-20.0,10.0,15.0,20.0]
print,'frequencia=', freq/1.e9

for npi=0,3 do begin
for iphi0=0,3 do begin
   phi0=phiv[iphi0]

   for iazim0=0,2 do begin
      azim0=azimv[iazim0]
      for ib_lat0=0,3 do begin
         b_lat0=b_latv[ib_lat0]
         r_arc0=0.007
            for ih=0,2 do begin
            par=h[ih]
         ofile='~/Dropbox/old/'+loop_file(st=st)
           
            freq=fmem
            nel=3.0e7
            n=32
            np=1.e8
            energy=[10.,1.e5]
            delta=[2.2]
            mmag=3000.0
            if file_search(ofile) eq '' then $
            calc_image,imagem,nel=nel,phi=phi0,b_lat=b_lat0,azim=azim0,freq=freq,$ 
               np=np,energy=energy,delta=delta, f_save=f_save, im=im, $ 
               mostra=mostra, r_sun=r_sun, mmag=mmag, par=par, $
               l_arc=l_arc, R_arc=R_arc,h_arc=h_arc, foot_s=foot_s, $
               px2cm=px2cm, mag_res=mag_res, ofile=ofile, ind=npi $ 
               ; print,ofile $
            else print,'Ja existe',' ',ofile
pula:
            endfor
       ;  endfor
      endfor
   endfor
endfor
endfor



end
