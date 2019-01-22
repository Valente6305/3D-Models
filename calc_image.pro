pro calc_image,imagem,st=st,nel=nel,phi=phi,b_lat=b_lat,azim=azim,freq=freq,$ 
               np=np,energy=energy,delta=delta, f_save=f_save, im=im, $ 
               mostra=mostra, r_sun=r_sun, mmag=mmag, par=par, $
               l_arc=l_arc, R_arc=R_arc,h_arc=h_arc, foot_s=foot_s, $
               px2cm=px2cm, mag_res=mag_res, ofile=ofile,nobrem=nobrem, $
               ind=ind, t=t, dir=dir,a_inc=a_inc,assim=assim

 ;; HISTORY:
 ;; JERCosta 01-march-2012 Created - version 1
 ;; PJSimoes 02-may-2012, added st_sun keywords to calc_image to
 ;; change dipole geometry
 ;; Add keyword for mag_res
 ;; JERCosta 13-may-2012, added keywords ofile, n, mmag

;Exemplo:  
;calc_image,imagem,phi=80.,azim=-5.0,b_lat=20,freq=[2.,3.75,9.4,17.0,35.]*1.e9
cd,curr=curdir
time0=systime(/sec)

if not keyword_set(dir) then dir=curdir
if not arg_present(st) then begin
 if not keyword_set(nel) then nel=3.0e7
 if not keyword_set(n) then n=16
 if not keyword_set(phi) then phi=30.
 if not keyword_set(b_lat) then b_lat=20.
 if not keyword_set(azim) then azim=20.
 if not keyword_set(freq) then freq=[3.75,9.4,17.0,35]*1e9
 if not keyword_set(np) then np=2.e8
 if not keyword_set(energy) then energy=[10.,5.e3]
 if not keyword_set(delta) then delta=[2.7]
 if not keyword_set(mag) then mag=3000.0
 if not keyword_set(T) then T=1.e8
 if not keyword_set(a_inc) then inc=0.0
 if not keyword_set(d_ang) then as=10.0
 if not keyword_set(par) then par=1.2 
 if not keyword_set(ind) then ind=0
 if keyword_set(nobrem) then bremss=0
endif else begin
   n=st.n
   nel=st.nel
   phi=st.phi
   b_lat=st.b_lat
   azim=st.az
   freq=st.freq 
   np=st.np ;1.e8 ; Dens. that will be changed by fill_np
   energy=st.en
   delta=st.delta
   im=st.m 
   mmag=st.B
   par=st.ie
   l_arc=1.0            ;Sempre UM
   R_arc=st.R_arc
   h_arc=st.h_arc
   foot_s=st.foot_s
   mag_res=mag_res
   ofile=loop_file(st=st)       
   nobrem=0             ;Sempre 0
   ind=st.ip       
   t=st.t
   inc=st.in
   as=st.assim
endelse

st_sun,st_arc,/new,show=0,azim=azim,b_lat=b_lat,/zoom,phi=phi,mag=mmag $
        , r_sun=r_sun,l_arc=l_arc, R_arc=R_arc,h_arc=h_arc, foot_s=foot_s, $
        a_inc=inc,assim=as,dens=np,nel=nel



;img=tvrd(/true)
; Para testes
;save,file='temp.sav',st_arc,n, np, t, mag_res, par, leng, nel

 mag_sun,st_arc,n, n_p=np, t=t, mag_res=mag_res, par=par, leng=leng, nel=nel ;,/debug

   ; Coisas para debug
 if keyword_set(debug) then begin
    pwp=where(mag_res.b ne 0)
    berror=(max(mag_res.b[pwp])-mag_res.st_arc.b*(1.+abs(mag_res.st_arc.assim)))/$
           mag_res.st_arc.b*(1.+abs(mag_res.st_arc.assim))
    print,'Bmax=',max(mag_res.b[pwp],min=mm),'   Bmin=',mm,'   Bmax(devido)=',$
          mag_res.st_arc.b*(1.+abs(mag_res.st_arc.assim)),'   Binput=',mag_res.st_arc.b,$
          '   Berror=',berror,'   Inc=',mag_res.st_arc.a_inc,mag_res.st_arc.d_ang
 endif
 print,'Repreenchimento do Np com fill_np: '
 print,'ip=',ind
 tt=fill_np(mag_res,ind)        ;preenche o ambiente com os swicthes
 if tt[0] ne -1 then mag_res.np=tt else  return 
; help,tt,np
; stop
; seed_file,st,np=max(mag_res.np)
;ofile=loop_file(st=st) 
 print,'Esse dir=>',dir
 print,'    file=>',ofile

transfer_2d_sol,delta, energy, mag_res.b, mag_res.theta, mag_res.nel, $
                 mag_res.np, freq, mag_res.px, imagem, transf, $
                 bremss=bremss, print=print, tau=tau, mag_res=mag_res, $
                 problem=problem,  f_save=f_save, temperature=t 

 px2cm=mag_res.px
 mag_res.st_arc.np=max(mag_res.np)
 t=where(mag_res.b ne 0)
 mag_res.st_arc.mirror=max(mag_res.b[t],min=mm)/mm

 if keyword_set(ofile) then begin
    flux=total(imagem.image[0].imagem)
    if flux lt 10. then save,file='~/Modelos_2016_Fluxo_baixos/'+ofile,delta,energy,par,ofile,mag_res,imagem  else save,file=dir+ofile,delta,energy,par,ofile,mag_res,imagem
    print,'Arquivo gravado:'    ;, ofile
 endif
 ftime=systime(/sec)-time0
 print,'Total time: ', ftime

 if keyword_set(mostra) then show_imagsp,imagem

end
 
