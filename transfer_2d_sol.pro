PRO transfer_2d_sol, delta, energy, bmag, angle, nel, np, freq, px, $
		     imagem, transf,  bremss=bremss, print=print, tau=tau, $
		     problem=problem, onestep=onestep, f_save=f_save, $
                     mag_res=mag_res, temperature=temperature
        ; JERCosta - Created - Jan-2012
	; Este transfer já calcula o gyro internamente e 
	; executa os passos da transferencia nos voxels
	; DELTA		- indice espectral da distribuicao dos eletrons em energia
	; ENERGY		- intervalo de energia em keV/MeV
	; BMAG		- e' um arranjo 3D das inducoes magneticas
	; ANGLE		- um arranjo simular a BMAG com os angulos de visada
	; NEL			- densidade numerica de eletrons em cm^-3
	; NP			- densidade do plasma em cm^-3
	; FREQ		- um escalar em Hz (alterado para vetor em 05/2010)
	; PX 			- e' a escala cm/pixel
	; IMAGEM		- e' o resultado imagem[3,*,*] => 0=I, 1=VX, 2=VO
	; BREMSS		- e' a keyword para incuir bremsstrahlung
	; TAU			- e' resultado
	; PROBLEM 	- e' 1 se deu problema e 0 se nao
	; ONESTEP 	- se igual a 1 faz transferencia em um passo (supoe homogeneo)

	; Constantes
	au=1.49597870D13				; unidade astronomica em cm
	dist=1.0*au
        dsource=px/au/!dtor*3600 ; 1.98	                ; precisamos verificar depois

	problem=0

	ss=size(bmag)
	IF ss[0] NE 3 THEN BEGIN			; testando se existe o arranjo 3D
	PRINT,'The magnetic field must be an array 3D'
		RETURN
	ENDIF
	sz=N_ELEMENTS(bmag)				; numero total de elementos no 3D

;	IF (SIZE(freq))[0] NE 0 THEN BEGIN		; testando se freq e' escalar
;		PRINT,'The freq must be one scalar only!'
;		RETURN
;	ENDIF

	; valores no centro do cubo para o caso homogeneo
	; p_B=B[ref.i,ref.j,ref.k]*1.d0                        
	; p_theta=theta[ref.i,ref.j,ref.k]*1.d0       

	imagem=DBLARR(2,ss[1],ss[2])			; definindo a imagem Io,Ix
	pp=WHERE(bmag NE 0,COMPLEMENT=pz)		; posicoes nao nulas do arranjo B
	zz=pp / (ss[1]*ss[2])				; posicoes z dos nao nulos (so' altura)
	mask=BYTE(bmag NE 0)				; 1 onde nao e' zero
	cc=TOTAL((bmag NE 0),3,/CUM)			; soma cumulativa dos 1's na mask
	IF (pz[0] NE -1) THEN cc[pz]=0 			; zero nos topos de col.s de materia
	bmag=REFORM(bmag,sz,/OVER)			; desmontando bmag
	angle=REFORM(angle,sz,/OVER)			; desmontando angle

	;st=size(T)					; nao existe temp por enquanto
	;if st[0] eq 0 then begin			; montando T
	;   T_mem=T
	;   T=mask*T
	;endif 

	; p_T = T[ref.i,ref.j,ref.k]			; nao existe temp por enquanto
	; T = REFORM(T,sz,/OVER)			; desmontando T

	sn = SIZE(nel)					; montando Nel
	IF sn[0] eq 0 THEN BEGIN
		nel_mem = nel
		nel = mask*nel 
	ENDIF
	; p_nel=Nel[ref.i,ref.j,ref.k]			; nao existe referencia
	nel = REFORM(nel,sz,/OVER)			; desmontando Nel

	sn=size(Np)
        np=REFORM(np,sz,/OVER)				; desmontando np
	;IF sn[0] eq 0 THEN BEGIN
	;   Np_mem=Np
	;   Np=mask*Np
	;endif else Np=mask*Np				; np e array
	;IF sn[0] EQ 0 THEN np_mem=np 			; np e escalar
	;np = mask*np
	;nnpp = REFORM(np,sz,/OVER)			; desmontando Np

	smax = 0.
	smax = MAX(cc)					; maximo percurso no volume
	pz = WHERE(cc EQ 0)				; onde nao existe fonte (materia)               
	s_image = TOTAL(mask,3)				; imagem do numero de pixels (2D)  
							; (s=depth para cada coluna)
	n_col = WHERE(s_image NE 0)			; indice das colunas nao nulas (2D)

	; ============== Geometria da transferencia================================
	; Determinacao no cubo de todas as colunas de emissao: 
	; 1) Onde estao os voxels emissores em cada nivel Z
	; 2) Qtas colunas existem em cada nivel ate' o nivel maximo smax
	; 3) Posicoes 3D colocadas em pos (XYZ). O count da' o n# de col. por nivel
	; 4) Posicoes 2D colcadas em ipos (XY naquele nivel)
	; =========================================================================
	pos = [0l]
	count = [0l]
	ipos = [0l]
	pppp = WHERE(s_image NE 0)
;print,'pppp=',pppp
	FOR i=0, smax-1l DO BEGIN
		ttt = WHERE(cc EQ (i+1), count0)
                sds=sort(ttt  mod  (ss[1]*ss[2]))
		IF i EQ 0 THEN tt0=MIN(ttt MOD (ss[1]*ss[2]),MAX=tt9)
                ttt=ttt[sds]
		pos = [pos,ttt]						; posicao em 3D (sequencial)
		hh = HISTOGRAM((ttt MOD (ss[1]*ss[2])),MIN=tt0,MAX=tt9)
		IF i EQ 0 THEN ppv = WHERE(hh NE 0)
		t_ipos = WHERE(hh[ppv] EQ 1,n_ipos) 
;print,'pos e ipos=',ttt mod 64,t_ipos
                ipos = [ipos,t_ipos]                            ; apenas a posicao dentro do arranjo 
                                                                ; de parameros enviados
		IF count0 NE n_ipos THEN BEGIN
		   PRINT, 'numero de colunas variou!'
		   STOP
		ENDIF
		count=[count,count0]					; qtos elementos por nivel
	ENDFOR
	count = count[1:*]
	inicio = TOTAL(count,/CUM)-count				; indice inicial de cada nivel
	fim = inicio+count-1						; indice final de cada nivel
	pos = pos[1:*]
	ipos = ipos[1:*]						; na origem do hh (sem tt0)
	level = pos MOD (ss[1]*ss[2])					; posicao XY no nivel

	IF KEYWORD_SET(onestep) THEN s_vec = s_image[n_col]*px ELSE s_vec=px

	; ============= T R A N S F E R E N C I A    R A D I A T I V A ==============
	; calculos para as bases das colunas
        ; Vamos calcular para todas as frequencias 
        nf=n_elements(freq)
        imt={imagem:fltarr(2,ss[2],ss[3])}
        image=replicate(imt,nf)
        in_freq=freq

        for i_freq=0, nf-1 do begin
           freq=in_freq[i_freq]
           indice0 = pos[inicio[0]:fim[0]]
           T1 = 0.
           T1 = SYSTIME(/SEC)
;stop
           M=0.
           gyro_trans_sol_im, delta, energy, nel[indice0], m, bmag[indice0]>0, np[indice0], $
                        angle[indice0], freq, res,bremss=bremss, px=px, problem=problem,temperature=temperature

           IF (problem EQ 1) THEN GOTO, erro
           transfer_kernel_sol, res, s_vec[0], transf, /reset, dsource=dsource, $
                                dist=dist, tau=tau ; COM RESET
;print,'EO nível           1=',res.res.eord, indice0 mod 64 ;;;transf.io
;stop
           transf_m=transf
;print,'tr nível           1=',transf_m.io,bmag[indice0],angle[indice0]
;print,pos[inicio[0]:fim[0]] mod 64
           indice0 = indice0 MOD (ss[1]*ss[2]) ; base da imagem
           im = DBLARR(ss[1],ss[2])

           IF KEYWORD_SET(onestep) THEN GOTO, pulo
           ;IF NOT KEYWORD_SET(onestep) THEN BEGIN
;print,transf_m[10].b,transf_m[50].b
;stop
           FOR i=1l,smax-1 DO BEGIN
              indice = pos[inicio[i]:fim[i]]
              transf = transf_m[ipos[inicio[i]:fim[i]]]
;stop
              gyro_trans_sol_im, delta, energy, nel[indice], m, bmag[indice]>0, np[indice], angle[indice], $
                              freq, rest,  bremss=bremss, px=px, problem=problem,temperature=temperature 
              IF (problem EQ 1) THEN GOTO, erro
              transfer_kernel_sol, rest, px, transf, dsource=dsource, $
                                   dist=dist, tau=tau_i ; SEM RESET
;print,'EO nível',i+1,'=',rest.res.eord,indice mod 64 ;ipos[inicio[i]:fim[i]]   ;transf.io
;stop
              transf_m[ipos[inicio[i]:fim[i]]] = transf
;print,'tr nível',i+1,'=',transf_m.io,bmag[indice],angle[indice]
;print,'pp nível',i+1,'=',pos[inicio[i]:fim[i]] mod 64, ipos[inicio[i]:fim[i]]
              tau[ipos[inicio[i]:fim[i]],*] = tau[ipos[inicio[i]:fim[i]],*]+tau_i
;if i eq 1 then stop
;print,transf_m[60].b,transf_m[50].b
           ENDFOR


           ;ENDIF

           ; =========================================================================
           pulo:

           im[pppp] = transf_m.Io*transf_m.o*transf_m.cgs2sfu
           imagem[0,*,*] = im
           im[pppp] = transf_m.Ix*transf_m.o*transf_m.cgs2sfu
           imagem[1,*,*] = im
                                ;im[pppp] = tau[*,0]
                                ;tau[*,0] = im
                                ;im[pppp] = tau[*,1]
                                ;tau[*,1] = im
           image[i_freq].imagem=imagem
           print,SYSTIME(/SEC)-t1
        endfor 
	bmag = REFORM(bmag,ss[1],ss[2],ss[3],/OVER)
	angle=reform(angle,ss[1],ss[2],ss[3],/OVER)
        np=reform(np,ss[1],ss[2],ss[3],/OVER)
	;IF sn[0] EQ 0 THEN nel=nel_mem ELSE
        nel=REFORM(nel,ss[1],ss[2],ss[3],/OVER)
        imagem={n:mag_res.n,st_arc:mag_res.st_arc,freq:in_freq,gyro:res.gyro, $
                lim:mag_res.lim,res:res.res,image:image}
	;STOP
	erro:
        if keyword_set(f_save) then save,file=f_save, imagem
	END

