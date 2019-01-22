pro grayplot,farr,xarr=xarr,yarr=yarr,string=string,sample=sample, $
	xtitle=xtitle,ytitle=ytitle

if(n_elements(xarr) eq 0 and n_elements(yarr) eq 0) then begin
	 nn=n_elements(farr(0,*))
   xarr=findgen(nn)
   yarr=findgen(nn)
endif

; create a dummy array for plotting
xsz=size(xarr)
ysz=size(yarr)
zarr=fltarr(xsz(1),ysz(1))*0.0

pltdev=!d.name
contour,[[0,0],[1,1]],/nodata, xstyle=4, ystyle = 4
px = !x.window * !d.x_vsize  ;Get size of window in device units
py = !y.window * !d.y_vsize 
swx = 1.*(px(1)-px(0))       ;Size in x in device units
swy = 1.*(py(1)-py(0))       ;Size in Y
aspi = 1.0       ;Image aspect ratio
aspw = swx / swy        ;Window aspect ratio
f = aspi / aspw         ;Ratio of aspect ratios
if f lt 1.0 then swx = aspi*swy
if f ge 1.0 then swy = swx/aspi
;swy=swx/xsz(1)*ysz(1)
;print,f,swx,swy

arr_sz=[swx,swy]
                       ; set dimensions
if pltdev eq 'PS' then begin
	arr_sz=256./arr_sz(0)*arr_sz
     if ((!p.multi(1)*!p.multi(2)) le 2.1) then arr_sz=512./arr_sz(0)*arr_sz
     print,'Gridding to array size ',arr_sz
;     !p.font=-1
endif

if(keyword_set(sample)) then $
tvscl,rebin(farr,xsz(1),ysz(1),/sample),px(0),py(0), $
             xsize=swx,ysize=swy,/device else $
tvscl,congrid(farr,arr_sz(0),arr_sz(1),interp=1,/center),px(0),py(0), $
             xsize=swx,ysize=swy,/device
contour,/nodata,zarr,xarr,yarr,/noerase, xstyle=1, ystyle = 1,$
        pos = [px(0),py(0), px(0)+swx,py(0)+swy],title=string,$
        /device,chars=0.8,ytitle=ytitle,xtitle=xtitle;,font=6
;xyouts,px(0)+swx/2,py(0)-swy/9.,string,/device,alignment=0.5

;if pltdev eq 'PS' then !p.font=-1

return
end
