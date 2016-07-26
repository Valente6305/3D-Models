function loop_file, st=st, file=file, dir=dir, freq=freq, new=new
;  Usage: # quero ver as modificações
;  name=loop_file(st=st) ===> simple return a name 
;  ttt=loop_file(file=file, st=st) --> then= st contain a structure 
ipar=0
ifile=0
inew=0
cd,current=curdir
if not keyword_set(dir) then dir=curdir
name=''
temst=0
temfile=0

   if arg_present(st) then temst=1
   if arg_present(file) then temfile=1
   if keyword_set(new) then begin
      st={n:32,phi:45.,B_lat:30.0,Az:15.0,in:0.0,assim:0.0,ie:0.34,ip:0,r_arc:0.02,$
          foot_s:0.02,h_arc:0.01,Nel:2.e7,np:1.e9,T:1.e8,En:[10.,100000], $
          delta:2.7,m:0,B:1000.,code:0,freq:[2.5,7.0]*1.e9}
   endif

   if keyword_set(freq) then begin 
      stt={n:st.n,phi:st.phi,B_lat:st.B_lat,Az:st.Az,in:st.in,assim:st.assim,ie:st.ie,ip:st.ip, $
           r_arc:st.r_arc,foot_s:st.foot_s,h_arc:st.h_arc,Nel:st.Nel,np:st.np, $
           T:st.T,En:st.En,delta:st.delta,m:st.m,B:st.B,code:st.code, $
           freq:freq}

      st=stt
   endif

if temst and temfile eq 0 then  begin
   ;v[0]=  n   (number of voxel)
   ;v[1]=  PHI (graus +/-)
   ;v[2]=  B_lat  (graus +/-)
   ;v[3]=  azim (graus +/-)
   ;v[4]=  in (graus +/-)
   ;v[5]=  par     (*100 elect. distr.)
   ;v[6]=  ip      (code for ambient dens. distr.)
   ;v[7]=  r_arc   (*1000. in solar radius)
   ;v[8]=  foot_s  (*1000. in solar radius) 
   ;v[9]=  h_arc   (*1000. in solar radius)
   ;v[10]=  alog10(nel)   (*100. and nel cm^-3)
   ;v[11]= alog10(temp)  (*100. and T in K)
   ;v[12]= energy[0]    (minimum in Kev)
   ;v[13]= energy[1]  (/1000. maximum in MeV)
   ;v[14]= delta     (*10. elct. spect. index)
   ;v[15]= m         (pitch angle distr. m=0 isotropic)
   ;v[16]= B    (maximum B in Gauss)
   ;v[17]= code  (TBD  = initially equal 0)
   if n_elements(st) eq 0 then begin   
      st={n:16,phi:45.,B_lat:30.0,Az:15.0,in:0.0,assim:0.0,ie:0.34,ip:0,r_arc:0.02,$
          foot_s:0.02,h_arc:0.01,Nel:2.e7,Np:1.e9,T:1.e8,En:[10.,100000], $
          delta:2.7,m:0,B:1000.,code:0, freq:[2.5,7.0]*1.e9}
   endif
;ost=st
   if n_tags(st) ne 20 then begin
      print,'Number of tags ne 20!'
      return, -1
      stop
   endif

   v=fltarr(20)
   v[0]=st.(0)
   v[1]=st.(1)
   v[2]=st.(2)
   v[3]=st.(3)
   v[4]=st.(4)
   v[5]=st.(5)*100.
   v[6]=st.(6)*100.
   v[7]=st.(7)
   v[8]=st.(8)*1000.
   v[9]=st.(9)*1000.
   v[10]=st.(10)*1000.
   v[11]=alog10(st.(11))*100.
   v[12]=alog10(st.(12))*10.
   v[13]=alog10(st.(13))*100.
   v[14]=(st.(14))[0]
   v[15]=(st.(14))[1]/1000.
   v[16]=st.(15)*10.
   v[17]=st.(16)
   v[18]=st.(17)
   v[19]=st.(18)
   tx=['nn','ph','la','az','in','as','ie','ip','c','fs','h','ln','lnp','lt',$
       'en','-','d','m','B','X','fr']
   fmt='(a2,i3,a2,i3,a2,i3,a2,i3,a2,i3,a2,i3,a2,i3,a2,i1,a1,i3,a2,i3,a1,i3,a2,i3,a3,i3,a2,' 
   fmt=fmt+'i3,a2,i3,a1,i3,a1,i3,a1,i1,a1,i4,a1,i1,a2)'

   novonome=string(tx[0],v[0],tx[1],abs(v[1]),tx[2],abs(v[2]),tx[3],abs(v[3]),$
                tx[4],v[4],tx[5],v[5],tx[6],v[6],tx[7],v[7],tx[8],v[8],tx[9],v[9], $
                tx[10],v[10],tx[11],v[11],tx[12],v[12],tx[13],v[13],$
                tx[14],v[14],tx[15],v[15],tx[16],v[16],tx[17],v[17],$
                tx[18],v[18],tx[19],v[19],tx[20],format=fmt)
   bw=byte(novonome)
   pp=where(bw eq 32)
   
   if pp[0] ne -1 then bw[pp]=48
   novonome=string(bw)

   if v[1] lt 0 then begin
      bw=byte(novonome)
      bw[7]=45
      novonome=string(bw)
   endif
   if v[2] lt 0 then begin
      bw=byte(novonome)
      bw[12]=45
      novonome=string(bw)
   endif
   if v[3] lt 0 then begin
      bw=byte(novonome)
     bw[17]=45
      novonome=string(bw)
   endif
   if v[4] lt 0 then begin
      bw=byte(novonome)
     bw[22]=45
      novonome=string(bw)
   endif

if v[5] lt 0 then begin
      bw=byte(novonome)
     bw[27]=45
      novonome=string(bw)
   endif

   nf=n_elements(st.freq)
   fmt='('+string(nf)+'i4)'
   w=string(st.freq/1e8,format=fmt)
   bw=byte(w)
   pp=where(bw eq 32)
   if pp[0] ne -1 then bw[pp]=48
   opar=novonome+string(bw)+'.sav'
   name=opar
   ipar=1
endif

;ost=st

if keyword_set(file) then begin
   pp=strpos(file,'nn')
   if pp[0] eq -1 then begin
      print, 'File does not begin with nn'
      return, -1
      stop
   endif
   pf=strpos(file,'.sav')
   pi=strpos(file,'fr')+2
   nf=(pf-pi)/4
   ost={n:32,phi:45.,B_lat:30.0,Az:15.0,in:0,assim:0,ie:0.34,ip:0,r_arc:0.02,foot_s:0.02,h_arc:0.01,Nel:2.e7,np:1.e9,T:1.e8,En:[10.,100000],delta:2.7,m:0,B:1000.,code:0, freq:fltarr(nf)}
   dir=strmid(file,0,pp)
   name=strmid(file,pp,255)
   ost.n=fix(strmid(file,pp+2,3))
   ost.phi=float(strmid(file,pp+7,3))
   ost.B_lat=float(strmid(file,pp+12,3))
   ost.Az=float(strmid(file,pp+17,3))
   if strmid(file,pp+22,2)  eq '--' then ost.in=float(-strmid(file,pp+24,1)) else ost.in=float(strmid(file,pp+22,3))
   if strmid(file,pp+27,2) eq '--' then ost.assim=float(-strmid(file,pp+29,1)) else  ost.assim=float(strmid(file,pp+27,3))
   ost.in=float(ost.in)
   ost.assim=float(ost.assim)
   ost.ie=float(strmid(file,pp+32,3))/100.
   ost.ip=fix(strmid(file,pp+37,1))
   ost.r_arc=float(strmid(file,pp+39,3))/1000.
   ost.foot_s=float(strmid(file,pp+44,3))/1000.
   ost.h_arc=float(strmid(file,pp+48,3))/1000.
   ost.Nel=10.^(float(strmid(file,pp+53,3))/100.)
   ost.np=10.^(float(strmid(file,pp+59,3))/10.)
   ost.T=10.^(float(strmid(file,pp+64,3))/100.)
   ost.En[0]=float(strmid(file,pp+69,3))
   ost.En[1]=float(strmid(file,pp+73,4))*1000.
   ost.delta=float(strmid(file,pp+77,3))/10.
   ost.m=fix(strmid(file,pp+81,1))
   ost.B=float(strmid(file,pp+83,4))
   ost.code=fix(strmid(file,pp+88,1))
   fmt='('+strcompress(string(nf),/remove_all)+'a4)'
   sw=strarr(nf)
   reads,strmid(file,pi,nf*4),sw,format=fmt
   ost.freq=float(sw)*1.e8
   ofile=strmid(file,pp,255)
   st=ost
   name=ofile
   ifile=1
endif
;if ipar eq 1 then name=opar
return, name
end


