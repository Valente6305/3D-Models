pro seed_file,st, new=new, nn=nn, phi=phi,b_lat=b_lat,in=in,assim=assim, $
              ie=ie,ip=ip,r_arc=r_arc,foot_s=foot_s, np=dens, $
              h_arc=h_arc,azim=azim,nel=nel,t=t,energy=energy, $
              delta=delta, mm=mm, mag=mag, code=code, freq=freq 
   ; n   (number of voxel)
   ; PHI (graus +/-)
   ; B_lat  (graus +/-)
   ; azim (graus +/-)
   ; in = inclinacao em graus 
   ; assim = assimetria do dipolo
   ; ie = par     (*100 elect. distr.)
   ; ip      (code for ambient dens. distr.)
   ; r_arc   (*1000. in solar radius)
   ; foot_s  (*1000. in solar radius) 
   ; h_arc   (*1000. in solar radius)
   ; alog10(nel)   (*100. and nel cm^-3)
   ; alog10(temp)  (*100. and T in K)
   ; energy[0]    (minimum in Kev)
   ; energy[1]  (/1000. maximum in MeV)
   ; delta     (*10. elct. spect. index)
   ; m         (pitch angle distr. m=0 isotropic)
   ; B    (maximum B in Gauss)
   ; code  (TBD  = initially equal 0)
;=================================================================
;                      Begin  Structure
;=================================================================
   if keyword_set(new) then begin
      tt=loop_file(st=st,freq=[3.75,9.4,17.0,35.0]*1.e9,/new)
   endif 
;=================================================================
;                      Change of Parameters
;=================================================================
   if size(nn,/type) ne 0 then st.n=fix(nn)
   if size(phi,/type) ne 0 then st.phi=float(phi)
   if size(b_lat,/type) ne 0 then st.b_lat=b_lat 
   if size(in,/type) ne 0 then st.in=in
   if size(assim,/type) ne 0 then st.assim=assim
   if size(ie,/type) ne 0 then st.ie=ie 
   if size(ip,/type) ne 0 then st.ip=ip
   if size(R_ARC,/type) ne 0 then begin 
      st.R_ARC=R_ARC 
   endif 
   if size(foot_s,/type) ne 0 then begin
      st.foot_s=foot_s
   endif
   if size(H_ARC,/type) ne 0 then begin 
      st.H_ARC=H_ARC 
   endif 
   if size(Azim,/type) ne 0 then st.az=azim
   if size(Nel,/type) ne 0 then begin 
      st.Nel=Nel
   endif
   if size(dens,/type) ne 0 then  st.np=dens ;acrecentei densidade 
      
   if size(T,/type) ne 0 then st.T=T
   if size(energy,/type) ne 0 then begin
      if n_elements(energy) ne 2 then begin
         print,'Energy need vector of two elements'
         print,'Energy not changed=>En=',st.en
      endif
      st.en=energy
   endif
   if size(delta,/type) ne 0 then st.delta=delta 
   if size(mm,/type) ne 0 then st.m=mm
   if size(mag,/type) ne 0 then begin 
      st.B=mag
   endif
   if size(code,/type) ne 0 then st.code=code
   if size(freq,/type) ne 0 then begin
      if n_elements(where(freq lt 1e4)) eq n_elements(freq) then freq=freq*1.e9
      stt={n:st.n,phi:st.phi,B_lat:st.B_lat,Az:st.Az,in:st.in,assim:st.assim,ie:st.ie,ip:st.ip, $
           r_arc:st.r_arc,foot_s:st.foot_s,h_arc:st.h_arc,Nel:st.Nel, np:st.np,$
           T:st.T,En:st.En,delta:st.delta,m:st.m,B:st.B,code:st.code, $
           freq:freq}

      st=stt
   endif
;=================================================================

end
