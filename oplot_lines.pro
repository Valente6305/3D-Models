;+
; PROJECT: w_mag
;
; NAME: oplot_lines
;
; PURPOSE: plot magnetogram selection with coordinates 
; 			  and overplot force free extrapolated lines
;			  as obtained from w_mag
;
; CATEGORY: Solar Physics
;
; CALLING SEQUENCE: oplot_lines,ST_LINES,ST_HEADER,RANGE[,lines=lines,/PS]
;
; INPUTS:
;				ST_LINES: Structure with magnetic field lines coordinates
;                     as generated by w_mag
;				ST_HEADER: Structure with magnetic information as generated
;                      by w_mag
;
; OPTIONAL INPUTS: 
;				LINES: Numeric array containing the index number of the lines
;						 to be overplotted. Equal to all if not set.
;
; KEYWORDS PARAMETERS:
;				PS: Set this keyword to send plots to a ps file named 'oplot.ps' at the session path
;
; OUTPUTS:   
;
; MODIFICATION HISTORY: J.E.R.Costa, 2010, v1.0
;		Keyword ps para gerar arquivo PS, Tereza Satiko, 2011
;		Overplot da seleção de linhas, Tereza Satiko, 2011
;-
;---------------------------------------------------------------

PRO oplot_lines,ST_LINES,ST_HEADER,RANGE,lines=lines,ps=ps

   soho_fac=0.99009902
	corner=st_lines.par.corner
   xc=st_lines.par.xc
   yc=st_lines.par.yc

	;loadct,33,file='magnetic.tbl',/silent
	loadct,33,file='~/idl/Magnetic/magnetic.tbl',/silent
	imagem=st_lines.par.img

   ;data=imagem
   si=SIZE(imagem)
   i_0=corner[0]
   i_9=corner[1]
   j_0=corner[2]
   j_9=corner[3]
   box_size=[i_9-i_0+1l,j_9-j_0+1l]
   box_center=[MEAN([i_0,i_9]),MEAN([j_0,j_9])]
   t_x=(FINDGEN(box_size[0])-(box_size[0]-1.)/2.+box_center[0])#REPLICATE(1,box_size[1])
   t_y=REPLICATE(1,box_size[0])#(findgen(box_size[1])-(box_size[1]-1.)/2.+box_center[1])
   tmp_img=INTERPOLATE(imagem,t_x,t_y,CUBIC=-0.5,MISS=0.0) ; zero fora da matriz
   ;era data[i_0:i_9,j_0:j_9]
   trange=range[[1,0,3,2]]

	IF NOT KEYWORD_SET(lines) THEN lines=LINDGEN(N_ELEMENTS(st_lines.nf)-1l)

	IF NOT KEYWORD_SET(ps) THEN BEGIN
		WINDOW,0,xsize=700,ysize=700
	   show_mag,tmp_img,xarr=trange[0:1], $
               yarr=trange[2:3], xtitle='<E-W> (arcsec)', ytitl=')cescra( vS-N^'
	   FOR i=0,N_ELEMENTS(lines)-1 DO BEGIN
			IF (lines[i] EQ 0l) THEN BEGIN
				inicio=0l 
				fim=st_lines.nf[0]-1l
			ENDIF ELSE BEGIN
				inicio=st_lines.nf[lines[i]-1l]
				fim=st_lines.nf[lines[i]]-1l
			ENDELSE
			OPLOT,(st_lines.xf[inicio:fim]-(st_header.x0))*st_header.px*soho_fac,$
		  		   (st_lines.yf[inicio:fim]-(st_header.y0))*st_header.px*soho_fac, $
		  		   ;thick=1,color=i mod 14 +241
					thick=2,color=250
	   ENDFOR
	ENDIF ELSE BEGIN
		SET_PLOT,'ps'
		DEVICE, FILE='oplot.ps', /COLOR, BITS=256
	   show_mag,tmp_img,xarr=trange[0:1], $
               yarr=trange[2:3],xtitle='<E-W> (arcsec)',ytitle='<N-S> (arcsec)'
	   FOR i=0,N_ELEMENTS(lines)-1 DO BEGIN
			IF (lines[0] EQ 0l) AND (i EQ 0l) THEN BEGIN
				inicio=0l
				fim=st_lines.nf[0]-1l
			ENDIF ELSE BEGIN
				inicio=st_lines.nf[lines[i]-1l]
				fim=st_lines.nf[lines[i]]-1l
			ENDELSE
			OPLOT,(st_lines.xf[inicio:fim]-(st_header.x0))*st_header.px*soho_fac,$
		  		   (st_lines.yf[inicio:fim]-(st_header.y0))*st_header.px*soho_fac, $
		  		   ;thick=2,color=i mod 14 +241
					thick=1,color=250
	   ENDFOR
		DEVICE, /CLOSE
		SET_PLOT, 'X'
	ENDELSE
END
