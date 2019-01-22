function mil_civtime, msecs, string=str, dsec=dsec, no_ms=no_ms

;;+
; NAME:
;	mil_civtime
; PURPOSE:
;	calculates civil time in hours, minutes, seconds, and
;	milliseconds from time in milliseconds
; CALLING SEQUENCE:
;	result = mil_civtime(msecs) 
; INPUTS:
;	msecs = any type except string and complex; contains time in
;	milliseconds of day.
; OPTIONAL INPUT PARAMETERS:
;	/string= if set time is returned as string : '14:46:23.456'
;	/no_ms = if set returned string does not contain ms : '14:46:23'
;	/dsec  = if set returned string contains only tenths of s :'14:46:23.5'
; OUTPUTS:
;	Function result = intarr(4) containing hours, minutes, seconds,
;	milliseconds unless keyword string is set.
; MODIFICATION HISTORY:
;	AM, Sep, 1990.
;	ER, Feb. 1996.  added keywords STRING, dsec
;-

 hour = fix(msecs/3.6d6)
 min  = fix(msecs/6d4 mod 60)
 sec  = fix(msecs/1d3 mod 60)
 msec = fix(msecs mod 1000)

if keyword_Set(str) then begin
   if keyword_set(no_ms) then begin
      sec=round(sec+msec/1000.)
      msec=0
   endif
   str=string(format='(2(i2.2,":"),i2.2)',transpose([[hour],[min],[sec]]))
   if max([msec ne 0]) then begin
      if keyword_set(dsec) then $
         str=str+string(format='(".",i1.1)',round(msec/100.)) $
      else str=str+string(format='(".",i3.3)',msec)
   endif
   return, str
endif
return, [hour,min,sec,msec]

end

; --------------------

pro mil_time_axis, milsecs, tickv, tickl, ticks, fine=fine, offset=offset

;;+
; NAME:
;        mil_time_axis
; PURPOSE:
;	 gives necessary information for plotting a time axis 
;	 labeled with civil time.
; CATEGORY:
;        plot utility
; CALLING SEQUENCE:
;        mil_time_axis, milsecs, tickv, tickl, ticks
; INPUTS:
;	 milsecs = 1 dim array of any type except string and complex; 
; 	 contains the times in milliseconds.
; OPTIONAL (KEYWORD) INPUT PARAMETERS:
;	 /fine = if set choose more ticks!
; OUTPUTS:
;        tickv = dblarr containing the values of the major ticks.
;	 tickl = string array containing labels for major ticks.
;	 ticks = lonarr(2) containing number of major and minor ticks.
; COMMON BLOCKS:
;        None.
; SIDE EFFECTS:
;	 None.
; RESTRICTIONS:
;	 None.
; MODIFICATION HISTORY:
;        AM, Sep, 1990. Modified by MG, Feb, 1991, ER, 1996.
;-

if n_elements(offset) eq 0 then offset=0l    ; IDL-BUG: plot only uses float!
tickv = make_array(30,/long  ,value=0)
tickl = make_array(30,/string,value=" ")
ticks = make_array(2,/long)

im   = n_elements(milsecs)-1
tb   = long(milsecs(0)+offset)  &  te = long(milsecs(im)+offset)
tlen = te-tb
if keyword_set(fine) then tlen=tlen/2 else fine=0

case 1 of 
  tlen ge 86400000L : begin & incr=21600000L& nt= 6 & nc= 1 & end  ;  6 h ticks
  tlen ge  7200000L : begin & incr=3600000L & nt= 6 & nc= 2 & end  ;  1 h ticks
  tlen ge  1200000L : begin & incr= 600000L & nt=10 & nc= 5 & end  ; 10 m ticks
  tlen gt   120000L : begin & incr=  60000L & nt= 6 & nc= 5 & end  ;  1 m ticks
  tlen gt    20000L : begin & incr=  10000L & nt=10 & nc= 8 & end  ; 10 s ticks
  tlen gt     2000L : begin & incr=   1000L & nt=10 & nc= 8 & end  ;  1 s ticks
  tlen ge      200L : begin & incr=    100L & nt=10 & nc=12 & end  ; .1 s ticks
  tlen ge       20L : begin & incr=     10L & nt=10 & nc=12 & end  ; 10ms ticks
  else :          message, 'time axis too short'
endcase

;
; construct labels
;
it    = tb/incr * incr
nlbl  = (te-tb)/incr +1 
nskip = nlbl/6 +1 - long(fine)
if nc eq 1 then nskip = (nskip - 1) >1

;

i = 0
while (it le te) do begin
   th = it/incr * incr 
   ctme = mil_civtime(th)
   neg=0
   if total(ctme) lt 0 then begin 
      ctme=abs(ctme)
      neg=1
   endif
   case nc of 
      1 :  str=string(format='(i0)',ctme(0) mod 24)   
      2 :  str=string(format='(i2.2)',ctme(0))   
      5 :  str=string(format='(i2.2,":",i2.2)',ctme(0), ctme(1))
      8 :  str=string(format='(2(i2.2,":"),i2.2)', ctme(0),ctme(1),ctme(2))
     12 :  str=string(format='(2(i2.2,":"),i2.2,",",i2.2)', $
                      ctme(0),ctme(1),ctme(2),fix(ctme(3)/10))
   endcase
   if neg then str='-'+str
   if (it ge tb) then begin
     tickv(i) = th - offset
     tickl(i) = str
     if (i mod nskip eq 0) then tickl(i) = str else tickl(i) = " "
     i=i+1
   endif
   it = it + incr
endwhile

ticks(0) = i-1
ticks(1) = nt
while (tickv(ticks(0)) gt te) do begin 
   ticks(0)=ticks(0)-1
   tickv(ticks(0)+1)=0.
endwhile



end

; --------------------------

pro t_plot, matrix,mmmm, error=error, $
          plot_error=plot_err,t_error=t_errorbar,ethick=errthick,eline=eline,$
          sm=smoothfaktor, ch=ch0, psym=psym, position=pos,$
          size=charsize, font=font, title=tit, $
          xmargin=xmargin,xtitle=xtit,xrange=xr,xrm=xrm, xa=xa, dx=dx,$
          ymargin=ymargin,ytitle=ytit,yrange=yr, ytype=ytype, ya=ya,$
          ytickv=ytickv,ytickn=ytickn,yticks=yticks, ytickformat=ytickf,$
          dy=dy,yscal=yscal,yminor=ymin, ystyle=ystyle, $
          noerase=noerase, overplot=overplot, color=color, thick=thick, $
          linestyle=linestyle, ynoclip=noclip, empty=empty, nsum=nsum, $
          r=r, g=g, b=b, fine=fine, min_value=min_value, max_value=max_value, $
          sst=sst, zero_offset=zero_offset,symsize=symsize,backg=backg,foreg=foreg,_EXTRA=ex

;+
; NAME:
;	t_plot
; PURPOSE:
;	plot time versus what ever you want
; CATEGORY:
;	plotting
; CALLING SEQUENCE:
;	t_plot, time, data     or    t_plot, lk
; INPUTS:
;	time = time array in ms since 00:00:00 UT
;	data = y-values
;	lk   = structure with first tag: time, following tags: data  OR
;		lk(0,*) = time,  lk(ch,*) = data
; OPTIONAL INPUT:
;	sst  = if set, time is in 0.1 ms instead of 1 ms, as used in SST data!
;	empty= if set, plot only the axis, no data!
;	error= array of the errorbars, if set only bars are plotted
;	sm   = smoothfaktor
;	ch   = channel: e.g. first data tag, second data tag, ... [ default=0 ]
;		if lk is a structure, ch can also be the tag name (string)
;	xrm  = xrange in minutes   [ the hour is taken from time(0) ]
;	/fine= if set the time axis has more ticks (finer scal)
;	xrange=xrange in ms since 0 UT or string: '13:15:37' or '13:14:51.572'
;	ya   = yticks:  1 = no ticks,           2 = only left, 
;			3 = only right,         4 = no ticks, large margins
;	xa   = x-axis:  1 = no x-axis, 	       -1 = x-axis no labels,
;			2 = bottom & labels, 	3 = bottom no labels, 
;			4 = top & labels, 	5 = top no labels.
;	plot_err= if set and err is given then also connect values.
;	ethick  = thickness of the errorbars [default=thick resp. !p.thick]
;	eline   = linestyle of errorbars [default=keyword linestyle, resp. 0]
;	t_error = if set plot double T errorbars (t_error is the width in %)
;	noerase = if set plot new data keeping the same x-axis with new y-axis.
;	overplot= if set plot new data on the same axis.
;	dx      = shift the x-axis by "dx" seconds
;	dy      = shift the plot by "dy" in the y-axis
;	yscal   = scale y-values:  data/yscal is plotted
;	ytype   = 0: linear y-axis,  1: log y-axis.
;	ynoclip = noclip keyword for plot, but only in y-direction.
;       foreg   = plot's foreground color (used on axis, titles and those
;                 infos that don't belong to the curve).
;       _EXTRA  = this keyword collects otherwise undefined keywords that
;                 are passed to the PLOT routine.
;
;	following keywords are forwarded to PLOT :
;	   psym, size=charsize, font,color, position, linestyle, thick, nsum
;	   title, xtitle, ytitle, xmargin, ymargin, ystyle
;	   min_value, max_value, yrange, ytype
;	   ytickv, ytickn, yticks, ytickformat, yminor, backg
;
; MODIFICATION HISTORY:
;	ER, IAP Uni Bern,  01/02/1995, added many keywords 1995/6
;	last change: 2 Oct 1996
;       Keyword improvements by J. GENTINA, CRAAM, October/2006
;-

on_error,2

if keyword_set(r) then begin & tvlct,255, 50,  0,1 & color=1 & endif
if keyword_set(g) then begin & tvlct,255,255,  0,2 & color=2 & endif
if keyword_set(b) then begin & tvlct,100,100,255,3 & color=3 & endif
if n_elements(ystyle)       eq 0 then ystyle=0
if not keyword_set(symsize) then symsize=1.0

ytl=!p.ticklen
if !y.ticklen ne 0 then ytl=!y.ticklen
if n_elements(ytit)         eq 0 then ytit= ''
yst=16 or ystyle
if n_elements(yr)           eq 0 then yr=[0,0] else yst=1 or ystyle
if keyword_set(ya)               then begin
      if ya gt 1 then if n_elements(xmargin) eq 0 then xmargin=[10,10]
      yst=(yst or 4)
      if keyword_set(ytickn) then ytn=ytickn else ytn=''
endif else ya=0
xst=1
if keyword_set(xa)               then xst=5 else xa=0
if xa eq -1                      then xst=1
if n_elements(pos)          eq 0 then pos=[0,0,0,0]
if n_elements(ch0)          eq 0 then ch0=0
if n_elements(dy)           eq 0 then dy=0
if n_elements(color)        eq 0 then color=(!d.n_colors-1)*(!D.name ne 'PS')
if n_elements(linestyle)    eq 0 then linestyle=0
if n_elements(eline)        eq 0 then eline=linestyle
if n_elements(ytickv)       eq 0 then ytickv=[0,0]
if n_elements(charsize)     eq 0 then charsize=1
if n_elements(xmargin)      eq 0 then xmargin=!x.margin
if n_elements(ymargin)      eq 0 then ymargin=!y.margin
if n_elements(ytype)        eq 0 then ytype=0
if n_elements(yticks)       eq 0 then yticks=0
if n_elements(ymin)         eq 0 then ymin=!y.minor
if n_elements(ytickn)       eq 0 then ytickn=''
if n_elements(ytickf)       eq 0 then ytickf=''
if n_elements(psym)         eq 0 then psym=0
if n_elements(thick)        eq 0 then thick=!p.thick
if n_elements(errthick)     eq 0 then errthick=thick
if n_elements(smoothfaktor) eq 0 then smoothfaktor=0

ch=ch0
sz=size(matrix)
if sz(sz(0) +1) eq 8             then begin
      szch = size(ch)
      if szch(szch(0)+1) eq 7 then begin
         strl=strlen(ch)
         ch = (where(strmid(tag_names(matrix),0,strl) eq strupcase(ch)))(0)-1
      endif
      chname=(tag_names(matrix))(ch+1)
      if ytit eq '1' then ytit=chname
      time=matrix.time 
      r   =matrix.(ch+1)
endif  else  begin
      if keyword_set(mmmm)    then  begin  
          time=reform(matrix)
          szm=size(mmmm)
          if szm(0) eq 2 and szm(1) gt 1 then r=reform(mmmm(ch,*)) else $
             r =reform(mmmm)
      endif else begin
          time=reform(matrix( 0,*))
          r   =reform(matrix(ch+1,*))
      endelse
endelse
if keyword_set(sst)   then time=time/10d
if keyword_set(yscal) then r=r/yscal
if keyword_set(error) then begin
   se=size(error)
   case se(0) of
      0:     err=matrix.err
      2:     err=reform(error(ch,*))
      else:  err=reform(error)
   endcase
   if smoothfaktor ge 2 then err=err/sqrt(smoothfaktor*1.)
   if keyword_set(yscal) then err=err/yscal
endif


if keyword_set(xrm)  then xr=long(time(0))/3600000L*3600000+xrm(0:1)*60000
if keyword_set(dx)               then time=time+dx*1000l
sxr=size(xr)
if sxr(sxr(0)+1) eq 7            then xr=time2ms(xr)
if n_elements(xr) eq 0           then xr=[min(time),max(time)]
if n_elements(xtit)         eq 0 then xtit= 'UT'
if n_elements(tit)          eq 0 then tit = ''
if smoothfaktor             ge 2 then r= smooth2(r,smoothfaktor)

p0=!p.position
!p.position=pos
if keyword_set(noerase) then begin
    !p.multi(0)=!p.multi(0)+1
    xr=!x.crange
endif

; IDL has a BUG that it transforms all arrays to float for plotting
; -----------------------------------------------------------------
common t_plot, offset
if not keyword_set(offset) then offset=0
if not keyword_set(overplot) and not keyword_set(noerase) then begin
   if keyword_set(zero_offset) then offset=0 else $
      offset = long(xr(0))/3600000L*3600000 
endif
time = time - offset
xr1 = xr
if not keyword_set(noerase) then  xr1 = xr - offset


if not keyword_set(overplot) or keyword_set(noerase) then begin
  mil_time_axis, xr1, tickv, tickl, ticks, fine=fine, offset=offset
     if xa mod 2 then begin tickl=strarr(29)+' ' & xtit='' & endif
     plot, time, r, xst=xst, yst=yst, xr=xr1, yr=yr, ytype=ytype, $
           xtickv= tickv,      xticks=ticks(0), $
           xtickname=tickl,    xminor=ticks(1), $
           xmargin=xmargin,    ymargin=ymargin, $
           xtitle= xtit, ytitle=ytit, title=tit,$
           yticks=yticks, ytickn=ytickn,ytickv=ytickv, $
           ytickf=ytickf, yticklen=ytl,yminor=ymin,$
           min_value=min_value, max_value=max_value, $
           charsize=charsize, font=font, backg=backg, col=foreg,/nodata,_EXTRA=ex
     if xa gt 1 then axis, xaxis=xa/4, /xst, xr=xr1, xtickn=tickl, $
             xtickv= tickv, xticks=ticks(0), xminor=ticks(1), xtit=xtit, $
             charsize=charsize, font=font,col=foreg
endif 

if ya gt 1 and ya ne 4 then  axis,yaxis=ya eq 3, /yst, yr=yr, ytickn=ytn, $
            yminor=ymin,ytickf=ytickf, yticks=yticks, ytickv=ytickv,ytit=ytit,$
            charsize=charsize, font=font
if ya gt 0 and ya ne 2 then plots,!x.crange(0),!y.crange,thick=!y.thick
if ya gt 0 and ya ne 3 then plots,!x.crange(1),!y.crange,thick=!y.thick
if keyword_set(empty) then return

if n_elements(noclip)       eq 0 then noclip=0 else begin
   if min(time) lt !x.crange(0) then begin
      r   = [interpol(r,time,!x.crange(0)),r(where(time gt !x.crange(0)))]
      if keyword_set(err) then  err = [0, err(where(time gt !x.crange(0)))]
      time= [!x.crange(0),time(where(time gt !x.crange(0)))]
   endif
   if max(time) gt !x.crange(1) then begin
      r   = [r(   where(time lt !x.crange(1))),interpol(r,time,!x.crange(1))]
      if keyword_set(err) then  err = [ err(where(time lt !x.crange(1))), 0 ]
      time= [time(where(time lt !x.crange(1))),!x.crange(1)]
   endif
endelse

if keyword_set(err) then begin
   for i=0l,n_elements(time)-1 do begin
      if err(i) ne 0 then plots,time(i),[-err(i),err(i)]+r(i)+dy,$
                color=color, linestyle=eline, thick=errthick, noclip=noclip
   endfor
   if keyword_set(t_errorbar) then begin
      wid=(!x.crange(1)-!x.crange(0))/200*[-t_errorbar,t_errorbar]
      for i=0,n_elements(time)-1 do begin
         if err(i) ne 0 then begin
            plots,time(i)+wid,-err(i)+r(i)+dy,col=color,lines=eline, $
                  thick=errthick, noclip=noclip
            plots,time(i)+wid, err(i)+r(i)+dy,col=color,lines=eline, $
                  thick=errthick, noclip=noclip
         endif
      endfor
   endif
endif 

if not keyword_set(err) or keyword_set(plot_err) then begin
   oplot,time,r+dy,psym=psym, col=color,line=linestyle,noc=noclip,$
         min_value=min_value, max_value=max_value, nsum=nsum,thick=thick,symsize=symsize
endif

!p.position=p0
end
