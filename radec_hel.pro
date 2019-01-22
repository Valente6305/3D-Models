;/*************************************************************
;*               RADEC_HEL                                    *
;*   RADEC_HEL conversts offsets in off_RA and off_DEC in     * 
;*   Heliographics  coordinnates Lat and Lon                  *
;*   off_ra and off_dec are in the same units as the radius R *
;*   i.e.: off_ra,off_dec and R are in arcmin                 *
;*         P, B, lat, lon are in deg
;**************************************************************
;*   J.E.R.Costa 09/10/06                                     *
;*                                                            *
;*                                                            *
;**************************************************************/

pro radec_hel, off_ra, off_dec, P,  B,  R, lat, lon

rad=!dtor ;

xlp=P*rad ;
xlp1=sin(xlp) ;
xlp2=cos(xlp) ;
xlb0=B*rad ;
xlb=-xlb0 ;
xlb1=sin(xlb) ;
xlb2=cos(xlb) ;
if ((R lt 15.72) and (R gt 16.7)) then R=16.4 ;
xlr=R ;
xlv1=off_dec/xlr ;  // arcmin/arcmin or other
xlv2=off_ra/xlr  ;  // arcmin/arcmin

xcom=xlv1*xlv1+xlv2*xlv2 ;
pp=where(xcom gt 1.0)
if pp[0] ne -1 then begin
  	print, 'There are coordinate(s) outside of solar disc' ;
	return ;
endif
xlv3=sqrt(1.00000000001-xcom) ;
xlxx=((xlv2*xlp1+xlv1*xlp2)*xlb2-xlb1*xlv3) ;
xlv4=atan(xlxx/sqrt(1.0-xlxx*xlxx)) ;
xlv5=cos(xlv4) ;
xlv6=xlv5*0.0
pp=where(xlv5 ne 0)
if pp[0] ne -1 then begin
	xlxx[pp]=((xlv2[pp]*xlp2[pp]-xlv1[pp]*xlp1[pp])/xlv5[pp]) ;
	xlv6[pp]=atan(xlxx[pp]/sqrt(1.0-xlxx[pp]*xlxx[pp])) ;
endif
lat=xlv4/rad   ;   //latitude South is negative
lon=-xlv6/rad  ;   //longitude East is negative

;// For the heliographic longitude the negative sign is East side ****

return ;

; // end of radec_hel
end
