pro roda_gyro
restore,'~/Dropbox/ProgramasIDL/Tratamento_NORH/status2.sav'
nomet='~/Dropbox/ProgramasIDL/Tratamento_NORH/25jan1/25jan1'
window,0,xsiz=960,ysize=540
for ii=233, 588-1 do begin
   indice=ii
   if ii eq 0 then indice=0.01
   if tab[0,ii] eq 2 then begin
      espectro, -1*indice,a, erro=er, oflux=oflux,ofreq=ofreq
      erro=er
      if erro then goto, pulaesse
      freq=[2.,3.75,9.4,17.0,35.]*1.e9   ;,80.]
      espec5,freq/1.e9,a,f
      dados=f     ;  espec5(freq,a)
      save,file='~/Dropbox/ProgramasIDL/Tratamento_NORH/meusdados.sav',dados,freq,a, oflux,ofreq
      print,'entrada: ',systime()
      pikaia_var, nomet, ii, result
      print,'saida: ',systime()
   pulaesse:
   endif
endfor
end
