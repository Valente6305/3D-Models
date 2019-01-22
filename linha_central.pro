pro linha_central,st_arc,f
f=abs(st_arc.assim)
;B0=32.*st_arc.r_sun^3*st_arc.deep^6*sqrt(4*st_arc.deep^2+st_arc.foot_s^2)*$
;     (st_arc.deep+st_arc.h_arc)^3/(sqrt(st_arc.deep^2+st_arc.foot_s^2)*$
;     (4*st_arc.deep^2+st_arc.foot_s^2)^(3))*st_arc.b
B0=st_arc.b0
window,0,xsiz=512,ysize=512
thetas=(findgen(200)/99.0-1.)*80.0*!dtor
desl=90.
xx=st_arc.r_sun*cos(thetas)-desl
yy=st_arc.r_sun*sin(thetas)
plot,xx,yy,xra=[0,40],yra=[-20.,20],/iso,yst=1 ; borda do Sol
r0=(st_arc.h_arc+st_arc.deep)*st_arc.r_sun
angf=atan(st_arc.foot_s/2,st_arc.deep) ; ang. polar da fotosfera
rf=r0*(cos(angf))^2                    ; raio vetor da fotosfera
ang=(findgen(400)/199.0-1.)*angf
r=r0*(cos(ang))^2
xx=r*cos(ang)+st_arc.r_sun-desl-st_arc.deep*st_arc.r_sun
yy=r*sin(ang)
oplot,xx,yy,thick=3   ; Loop central sem assimetria
xyouts,20,18,'r0='+string(r0),chars=1.3
xyouts,20,16,'rf='+string(rf),chars=1.3
xyouts,20,14,'Lf='+string(angf/!dtor),chars=1.3
xyouts,20,12,'B0='+string(b0),chars=1.3
xyouts,20,10,'Foot_S='+string(st_arc.foot_s),chars=1.3
xyouts,20,08,'Height='+string(st_arc.h_arc),chars=1.3
bf=b0*sqrt(1.0+3.*sin(angf)^2)/rf^3; campo mag nos pes
bt=b0/r0^3 ; Campo mag no topo
xyouts,20,06,'Bfoot='+string(bf),chars=1.3
xyouts,20,04,'Bfot(input)='+string(st_arc.b),chars=1.3
xyouts,20,02,'Btopo='+string(bt),chars=1.3
xyouts,20,0,'Mirror='+string(bf/bt),chars=1.3

      if st_arc.assim lt 0 then signal_assim=-1.0 else signal_assim=1.0
      k=((1./(1.+st_arc.assim))*rf^3/sqrt(1+3.*sin(angf)^2))^2
      r2=(findgen(301)/300.+0.5)*rf
      v=r2^6-k*(1.+3.d0*(1.-r2/r0))
      if (where(v*v[1:*] le 0))[0] eq -1 then stop
      tmp=min(abs(v),ii)
      r2=r2[ii]
      ang2=acos(sqrt(r2/r0))
      Bn=st_arc.b0*sqrt(1.0+3.*sin(ang2)^2)/r2^3
      a2=r2^2+(st_arc.deep*st_arc.r_sun)^2-$
         2.*r2*st_arc.deep*st_arc.r_sun*cos(ang2)
      phi=!pi/2.-asin(r2*sin(ang2)/sqrt(a2))
      phi=signal_assim*phi/!dtor
print,k,r2,bn,a2,ang2,phi



k=((1./(1.+f))*rf^3/sqrt(1+3.*sin(angf)^2))^2
r2=(findgen(301)/300.+0.5)*rf
v=r2^6-k*(1.+3.d0*(1.-r2/r0))
if (where(v*v[1:*] le 0))[0] eq -1 then stop
tmp=min(abs(v),ii)
r2=r2[ii]
ang2=acos(sqrt(r2/r0))

Bn=b0*sqrt(1.0+3.*sin(ang2)^2)/r2^3
xyouts,20,-2,'Bnovo='+string(bn),chars=1.3
xyouts,20,-4,'Fracao='+string((bn/st_arc.b-1.0)*100.),chars=1.3

a2=r2^2+(st_arc.deep*st_arc.r_sun)^2-2.*r2*st_arc.deep*st_arc.r_sun*cos(ang2)
phi=(!pi/2.-asin(r2*sin(ang2)/sqrt(a2)))/!dtor/!dtor

print,k,r2,bn,a2,ang2,phi

; st_arc.d_ang=phi/!dtor  !!!!!!!

xyouts,20,-6,'D_ang='+string(phi/!dtor),chars=1.3
T3D,/reset,translate=[-(st_arc.r_sun-desl),0.,0], $
    rotate=[0.,0.,phi/!dtor],matrix=tr0 
T3D,tr0,translate=[(st_arc.r_sun-desl),0.,0],matrix=tr0
c=[[[xx]],[[yy]],[replicate(0.,400)],[replicate(1.,400)]]
c=c#tr0
oplot,c[*,0],c[*,1],col=150

end
