pro resolve_flux

gyro,[5.,10.,20.,30.]*1.e9,fluxg,/plot,size=1,heig=4.e7,e1d=e1dg,e2d=e2dg,a1d=a1dg,a2d=a2dg,delta=3.0


gyrosyncs,3.,[10.,5e3],1.e7,0,500.,1.e9,45.,[5.,10.,20.,30.]*1.e9,fluxs,theta=1, $ 
          lsource=4.e7,e1d=e1ds,e2d=e2ds,a1d=a1ds,a2d=a2ds,/transf,/sfu

t1=total((e1dg-e1ds)/e1ds)
t2=total((e2dg-e2ds)/e2ds)
t3=total((a1dg-a1ds)/a1ds)
t4=total((a2dg-a2ds)/a2ds)
print,'Total emis o,x=' , t1,t2
print,'Total abs  o,x=' , t3,t4


stop
end
