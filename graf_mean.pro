pro graf_mean,x,y,bin,limite=limite,ytit=ytit,xtit=xtit
; Plot de um y vs x onde x é acumulado em bins (hostograma) 
; e y é uma média dentro de cada bin.

if not keyword_set(xtit) then xtit=''
if not keyword_set(ytit) then ytit=''
if not keyword_set(bin) then bin=1.
if not keyword_set(limite) then $
  hist = HISTOGRAM(x,omin=dmin,omax=dmax,bin=bin,reverse_ind=r) $
else $
  hist = HISTOGRAM(x,omin=dmin,omax=dmax,bin=bin, $
         min=limite[0],max=limite[1],reverse_ind=r) 
bins = FINDGEN(N_ELEMENTS(hist))*bin + dmin

yy=fltarr(n_elements(bins))
j=0l
for i=0, n_elements(bins)-1 do begin
   IF R[i] NE R[i+1] THEN begin 
      yy[j]=mean(y[r[R[I] : R[i+1]-1]])
;      print,'Bin ',j, x[r[R[I] : R[i+1]-1]]
;      print,'media=',yy[j]
;      print,'----------'
   endif else begin
      yy[j]=0.
   endelse
   j=j+1
endfor
PLOT, bins, yy, PSYM = 10, $
   XTITLE = xtit, YTITLE = ytit ,xran=[dmin,dmax],chars=1.5

end
