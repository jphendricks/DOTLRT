;plot sib3 soil characteristics
; such as
; thermal conductivity as a function of temperature 
; and moisture, for a given clay fraction, sand fraction,
; and porosity.
 
 
PRO sib3_soil_props, ps=ps
 
device, decomposed = 0
 
if keyword_set(ps) THEN BEGIN
 
   set_plot,'ps'
   !P.Font=0
 
   DEVICE,Bits_Per_Pixel=8, /Color, $
     Filename='CCSM_monthly_1.ps', $
    /Helvetica, /Bold, Font_size=9,          $
    XOffset=0.05, YOffset=0.5, XSize=7.90, YSize=9.5, /Inches
 
endif

 
;mess with some color table stuff
 
 
;color #  0  1  2  3  4  5  6  7   8  9 10  11  12
RED    = [0, 1, 1, 0, 0, 1, 0, .5, 0, 0, 1, .5, .9]
GREEN  = [0, 1, 0, 1, 0, 1, 1, 1, 0, .5, 0, .1, .9]
BLUE   = [0, 1, 0, 0, 1, 0, 1, 0, .5, 0, 1, .1, .9]
TVLCT, 255*RED,255*GREEN,255*BLUE
 
 
 
;display window stuff
if NOT keyword_set(ps) Then Window, XSize=620, YSize=800
if keyword_set(ps) Then plot_thick = 10 ELSE plot_thick = 4
 
!P.color=0
!P.background=1
!P.MULTI = [0,1,2,0,1]

;----HERE ARE THE VARIABLES YOU MODIFY (CURRENT VALUES FROM WLEF)

poros = 0.4424
cfrac = 0.15
sfrac = 0.37

fl = 1.0   ; FRACTION LIQUID FOR TEMPS<273

;----END MODIFICATION

; some variables
;-------------------

satw = FINDGEN(100)
satw(*) = satw(*) * 0.01 + 0.01

temp = FINDGEN(100)
temp(*) = temp(*) + 224.0

fliq = FINDGEN(100)
fliq(*) = fliq(*) *0.01

conductivity = DBLARR(100,100)

;-------------------


;---some vars dependent on soil chars (poros, sfrac, cfrac)

bd = (1.0 - poros) * 2700.0

tkm = (8.80 * sfrac + 2.92 * cfrac)/(sfrac + cfrac)

tkmg = tkm * (1.0 - poros)

tksatu = tkmg * 0.57^poros

tkdry = (0.135 * bd +64.7) / (2700.0 - 0.947 * bd)


;---loop over temperature and soil moisture fraction

for i = 0,99 do begin
  for j = 0,99 do begin

     if(temp(i) GT 273.16) then begin

;        print,i,j,satw(j),alog10(satw(j))+1.0
        dke = alog10(satw(j)) + 1.0

        if(dke LT 0.0) then dke = 0.0

        dksat = tksatu

     endif else begin
 
        dke = satw(j)

        dksat = tkmg * 0.249^(fl * poros) * 2.29^poros

    endelse

;     print,i,j,temp(i),satw(j),dke,dksat

     conductivity(i,j) = dke * dksat + (1.0 - dke) * tkdry

  endfor
endfor

;contour it up!

surface,conductivity

;contour, conductivity

for i = 0,99 do begin
  for j = 0,99 do begin


        dke = satw(j)

        dksat = tkmg * 0.249^(fliq(i) * poros) * 2.29^poros


     conductivity(i,j) = dke * dksat + (1.0 - dke) * tkdry

  endfor
endfor


surface, conductivity


if keyword_set(ps) THEN BEGIN
   Device,/CLOSE
   Set_Plot, 'X'
   !P.Font=-1
endif
 
 
 
 
end
