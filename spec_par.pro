function spec_par,s,area=area
; s= fluxo nas quatro frequencias
; x=[findgen(1000)-300.]/30.  ; Peso para alfa fino (x=-alfa_fino)
; plot,x,(x>0)/(.7+(x>0))     ; proximo de zero e positivo tem peso=0.
; Delta alfa_fino determina NL
; Cruzamento da reta espessa com fina da a freq_turnover
; Delta alfa_espesso da anisotropia
; spar=[alfa_espesso,turnover,alfa_fino,S34]
if not keyword_set(area) then tarea=[1.0,1.0] else tarea=area
ss=size(s)
freq=[3.75,9.4,17.0,35.0]  ; GHz
if ss[1] ne 4 then stop
b=alog10(s[1]/s[0])/alog10(freq[1]/freq[0])
d=alog10(s[3]/s[2])/alog10(freq[3]/freq[2])
a=alog10(s[1])-b*alog10(freq[1])
c=alog10(s[2])-d*alog10(freq[2])
alfa_esp=b
alfa_fin=-d
turn=(a-c)/(d-b)
s17=s[2]/tarea[0]*(alfa_fin > 0.)/(0.7+(alfa_fin > 0.))
s34=s[3]/tarea[1]*(alfa_fin > 0.)/(0.7+(alfa_fin > 0.))
if abs((atan(b)-atan(d))/!dtor) lt 1.0 then begin
   turn=mean(freq[1:2])
   if b ge 0 and d ge 0 then turn=freq[3] 
   if b lt 0 and d ge 0 then turn=freq[3] 
endif
if b lt 0 and d lt 0 then turn=freq[0] 
spar=[alfa_esp,turn,alfa_fin,s17,s34]

return,spar
end
