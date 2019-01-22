pro roda_loop

dir='~/Modelos/'
restore,'~/Dropbox/ProgramasIDL/Tratamento_NORH/theta.sav'

phi=80.0
b_lat=22.0

azimv=[-80.0,-55.,-30.,0.,30.,55.,80]

freq=[3.75,9.4,17.0,35.]*1.e9  ;
fmem=freq
nn=32
r_arc=0.007
foot_s=0.02
h_arc=0.02
np=1.e8
energy=[10.,1.e5]
code=0
mm=0
T=1.e8
h=[0.05,0.34,2.0]  ; preenchimento


print,'frequencia=', freq/1.e9

temporario=1
for ih=0,2 do begin
for i=0,6 do begin
; i=randomu(seed)*(n-1l)
 mag=1366.0
 delta=2.1
 nel=1.95e7
 phi0=phi
 b_lat0=b_lat
 npi=0   ; NAO MEXA AQUI 
      azim0=azimv[i]
      r_arc0=r_arc
;      ih=1            ; Aqui o preenchimento (0.05,0.34,2.0)
         ie0=h[ih]
         mostra=0
         seed_file,st, /new, nn=nn, phi=phi0,b_lat=b_lat0, $
              ie=ie0,ip=npi,r_arc=r_arc0,foot_s=foot_s, $
              h_arc=h_arc,azim=azim0,nel=nel,t=t,energy=energy, $
              delta=delta, mm=mm, mag=mag, code=code, freq=freq
         ofile=loop_file(st=st)
         freq=fmem
         if file_search(dir+ofile) eq '' then begin
            calc_image,imagem,st=st, f_save=f_save, mostra=mostra, $
               px2cm=px2cm, mag_res=mag_res, ofile=ofile, ind=npi, dir=dir
            print,ofile 
         endif else print,'Ja existe',' ',ofile
pula:
     



endfor
endfor
end
