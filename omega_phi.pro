;Bx = Magnetic field x component. 2D array.
;By = Magnetic field y component. 2D array.
;Bz = Magnetic field z component. 2D array. 
;B0 = sun center latitude in deg (B0 positive we see north pole!)
;     This definition +/- is different of the old definition!
;tr = Transform matrix for Bx,By,Bz in the viewer direction
;============================================================

function omega_phi,Bx,By,Bz,r_s,tr,coord,im,ang_v

s=size(im)
pp_bneg=where(im lt 0,complement=pp_bpos)
ma=max(im,min=mi)
mi=-mi

theta=ang_v
if pp_bneg[0] ne -1 then theta[pp_bneg]=180.-ang_v[pp_bneg]
im0=im/cos(ang_v*!dtor)
n=s[1]*s[2]

; phi angle is the angle between B(extrap.) and radial direction
b=sqrt(bx^2+by^2+bz^2)
cos_phi=bz/(b>0.01) < 1.0 > (-1.0)

; Omega is the angle between B and the viewer
coord0=[[coord],[replicate(1.,n)]]
b_g=(coord0#invert(tr)+[[reform(bx/b,n)],[reform(by/b,n)], $
    [reform(bz/b,n)],[replicate(0.,n)]])#tr-coord0
;coord=[[reform(bx/b,n)+r_s.x0],[reform(by/b,n)+r_s.y0],[reform(bz/b,n)+511.5],[replicate(1.,n)]]
;coord=coord#tr
;cos_omega=reform(coord[*,2]-511.5,s[1],s[2]) < 1.0 > (-1.0)
cos_omega=reform(b_g[*,2],s[1],s[2]) < 1.0 > (-1.0)

;Limitando o cos_omega se o im e' pequeno. Somente B mais ou menos 
; alinhado com Bl sao verdadeiros.
cos_omega[pp_bpos]=cos_omega[pp_bpos] > 0.7
cos_omega[pp_bneg]=cos_omega[pp_bneg] < (-0.7)

pp=where(im*cos_omega gt 0, complement=ppn)

if pp[0] ne -1 then begin
   factor=1.0/abs(cos_omega[pp])*(abs(cos_phi[pp])>0.7) ; limitado em phi tbem
   im0[pp]=im[pp]*factor
endif
print,n_elements(pp),'  in',n

return,im0

end
