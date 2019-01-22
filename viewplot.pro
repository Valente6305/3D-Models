pro viewplot
restore,'25jan1/25jan1_total.sav'
histog_j,tab_str.par[0],5,'Theta (deg)',res=res
;histog_j,tab_str.par[1],1,'Delta ',res=res
;histog_j,tab_str.par[2],0.1,'Nel (cm-3)',res=res
;histog_j,tab_str.par[3],0.1,'L (cm)',res=res
;histog_j,tab_str.par[4],25,'Bmag (Gauss)',res=res
;histog_j,tab_str.par[5],1,limite=[1.,60.],'Size (arcsec)',res=res
;Use: IF R[i] NE R[i+1] THEN A[R[R[I] : R[i+1]-1]] = xxxxx
r=res.r 
xx=''
;plot,findgen(1920)-960.0,findgen(1920)-960.0,/nodata 
plot,findgen(90),findgen(90)*10.,/nodata  
conta=1                                    
for i=0,n_elements(r)-3 do begin
   IF R[i] NE R[i+1] THEN begin 
;      plot,findgen(1920)-960.0,findgen(1920)-960.0,/nodata   
;      xyouts,tab_str[R[R[I] : R[i+1]-1]].x,tab_str[R[R[I] : R[i+1]-1]].y, $
      xyouts,tab_str[R[R[I] : R[i+1]-1]].thetac/!dtor,tab_str[R[R[I] : R[i+1]-1]].y, $
             strcompress(string(conta)),alig=0.5
      conta=conta+1
   endif
read,xx
endfor

end
