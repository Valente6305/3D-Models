pro file_distrib, dir, st, names=names
; Input: dir (where models are) 
; Outputs: st (structure with all propertie of models in Dir)
;          names (names of models in Dir)
dir='~/Modelos_ass'
if not keyword_set(names) then ff=file_search(dir +'/nn*.sav') else ff=names
names=ff
nf=n_elements(ff)

for i=0l, nf-1 do begin
   ttt=loop_file(file=ff[i],st=t)
   if i eq 0 then begin 
      tt={n:t.n,phi:t.phi,B_lat:t.b_lat,Az:t.az,in:t.in,as:t.as,ie:t.ie,ip:t.ip,r_arc:t.r_arc,$
          foot_s:t.foot_s,h_arc:t.h_arc,Nel:t.Nel,T:t.t,En:t.en, $
          delta:t.delta,m:t.m,B:t.b,code:t.code, freq:''}

      st=replicate(tt,nf)
      pp=strpos(ff[0],'fr')+2
      np=strpos(ff[0],'.sav')-1-pp+1
      st[0].freq=strmid(ff[0],pp,np)
   endif else begin
      pp=strpos(ff[i],'fr')+2
      np=strpos(ff[i],'.sav')-1-pp+1
      st[i]={n:t.n,phi:t.phi,B_lat:t.b_lat,Az:t.az,in:t.in,as:t.as,ie:t.ie,ip:t.ip,$
          r_arc:t.r_arc,$
          foot_s:t.foot_s,h_arc:t.h_arc,Nel:t.Nel,T:t.t,En:t.en, $
          delta:t.delta,m:t.m,B:t.b,code:t.code, freq:''}
      st[i].freq=strmid(ff[i],pp,np)
   endelse
endfor
;stop
end
