Function quadrante,m, rad=rad
if keyword_set(rad) then sc=!pi else sc=180.
pp=where(m gt sc/2.)
p=where(m eq 0)
n=m
if pp[0] ne -1 then n[pp]=sc-n[pp]
n[p]=0
return,reform(n)
end
