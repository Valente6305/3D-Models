pro ofile

ff=file_search('~/kawate/*.sav') 
for i=0l,n_elements(ff)-1 do begin
   restore,ff[i]
   r=strpos(ofile,'/',/reverse_searc)     
   w=strmid(ofile,r+1,255)                
   ;ofile=strmid(w,0,255)
   ;nome=strmid(w,0,25)+'as000'+strmid(w,25,255)
      print,ofile
      wait,0.0001
endfor
end
