pro plotcool

; Generates a 3-phase cooling curve and determines
; the equilibrium temperature-density relation (also, pressure-density),
; stored in neqi and Teqi (and peqi).
; An effective gamma ("geff") is also determined via d ln(P)/d ln rho.

  common wurstbrot,Tmin,Tmax,nmin,nmax,nshift,eps,p1shift

  nres   = 512l
  nmin   = 1d-5
  nmax   = 1d3
  Tmin   = 1d1
  Tmax   = 1d7
  p1shift= 1d0   ; location of 1st peak. For 1d0:npeak=1.5d-4; 0.3d0:npeak=4.8d-5,3d0:npeak=3.6d-4
  eps    = 1d-6
  nshift = 3d0 ; amount by which peak of P-n curve is shifted
   
;===================================================

  dens = 1d1^((alog10(nmax)-alog10(nmin))*dindgen(nres)/double(nres)+alog10(nmin))
  temp = 1d1^((alog10(Tmax)-alog10(Tmin))*dindgen(nres)/double(nres)+alog10(Tmin))
  neqi = dens          
  teqi = dblarr(nres) ; equilibrium temperature
  peqi = dblarr(nres) ; equilibrium pressure
  fcool= dblarr(nres)

  file = '/hpc/scratch/astro/users/jmg2223/proteus/run/hvc/idl/plotcool.dat'
; now find the equilibrium temperature for each density (simple bisearch)
;  window,0,xsize=600,ysize=400
;  plot,alog10([nmin,nmax]),alog10([Tmin,Tmax]),/nodata,xtitle='log n',ytitle='log T'
  for i=0,nres-1 do begin
    Tmid = getteqi(neqi[i])
;    printf,file,format='("PLOTCOOL: i = ",i4,", n = ",e13.5,", T = ",e13.5,", P = ",e13.5)',i,neqi[i],Tmid,neqi[i]*Tmid
    print,'(PLOTCOOL: i = ",i4,", n = ",e13.5,", T = ",e13.5,", P =,e13.5)',i,neqi[i],Tmid,neqi[i]*Tmid
    teqi[i] = Tmid 
    peqi[i] = Tmid*neqi[i]
;    oplot,alog10([neqi[i],neqi[i]]),alog10([teqi[i],teqi[i]]),psym=4,symsize=0.5
  endfor

  minp = 0.0;min(alog10(peqi))
  maxp = 5.0;max(alog10(peqi))

; effective gamma: dln P/dln rho
  geff    = dblarr(nres)
  for i=0,nres-1 do begin
    geff[i] = getgamma(neqi[i])
  endfor
  
  dgamdrho = dblarr(nres)
  for i=0,nres-1 do begin
    dgamdrho[i] = dgammadrho(neqi[i])
  endfor

;  window,1,xsize=600,ysize=400
;  plot,alog10(neqi),alog10(peqi),yrange=[minp,maxp],psym=4,symsize=0.5,xtitle='log n',ytitle='log P'
;  oplot,[-2.0,4.0],[-1.0,5.0],line=2  ; T=10K
;  oplot,[-5.0,1.0],[-1.0,5.0],line=2 ; T= 10^4K
;  oplot,[-7.0,0.0],[-1.0,6.0],line=2 ; T= 10^6K

;  window,2,xsize=600,ysize=400
;  plot,alog10(neqi),geff,psym=4,symsize=0.5,xtitle='log n',ytitle='gamma'

;  window,3,xsize=600,ysize=400
;  plot,alog10(neqi),dgamdrho,psym=4,symsize=0.5,xtitle='log n',ytitle='d gam d rho'

;===================================================

  return
end

;===================================================

function getcoolfunc,n,T
; see Inoue, Inutsuka & Koyama 2007, ApJL, 658, 99
  common wurstbrot
  mix1     =   0.5d0*(1d0+tanh((n-4.5d-3)/5.2d-4))
  mix2     =   0.5d0*(1d0-tanh((n-4.5d-3)/5.2d-4))
  if (mix1+mix2 ne 1.0) then print,mix1+mix2
  heatcool =   mix1*(   2d-26 $ ; a constant heating term 
                      - (n*nshift)*(2d-19*exp(-1.184d5/(T+1d3)) + 2.8d-28*sqrt(T)*exp(-9.2d1/T))) $; classical TI
             + mix2*(   p1shift*3.0d-29 $
                      - (n*nshift*20.0)*(2d-22*exp(-1.2d7/(T+1d5)) + 3.5d-30*sqrt(T)*exp(-(p1shift)^0.5*9.2d4/T))) 
             
  return,heatcool
end

;===================================================

function getcool,n,T
  common wurstbrot
  mix1     =   0.5d0*(1d0+tanh((n-4.5d-3)/5.2d-4))
  mix2     =   0.5d0*(1d0-tanh((n-4.5d-3)/5.2d-4))
  c0 =   mix1*(2d-19*exp(-1.184d5/(T+1d3)) + 2.8d-28*sqrt(T)*exp(-9.2d1/T)) $; classical TI
       + mix2*(2d-22*exp(-1.2d7/(T+1d5)) + 3.5d-30*sqrt(T)*exp(-9.2d4/T))
  return,c0
end

;===================================================
function getteqi,rho
  common wurstbrot
  Tlo  = Tmin
  Thi  = Tmax
  Tmid = 0.5d0*(Thi+Tlo)
  Pmid = 0d0
  ic   = 0l
  while (abs(Thi-Tlo)/Tmid gt eps) do begin
    c0 = getcoolfunc(rho,Tmid)
    if (c0 gt 0d0) then begin
      Tlo = Tmid
    endif else begin
      Thi = Tmid
    endelse
    Tmid = 0.5d0*(Thi+Tlo)
    ic = ic+1l
  endwhile
  return,Tmid
end

;===================================================
function getgamma,rho
  common wurstbrot
  epsrho = 1d-2
  teqip = getteqi(rho*(1d0+epsrho))
  teqim = getteqi(rho*(1d0-epsrho))
  geff  = 1d0+(alog(teqip)-alog(teqim))/(alog(rho*(1d0+epsrho))-alog(rho*(1d0-epsrho)))
  return,geff
end

;===================================================
function dgammadrho,rho
  common wurstbrot
  epsrho = 1d-2
  rhop  = rho*(1d0+epsrho)
  rhom  = rho*(1d0-epsrho)
  teqip = getteqi(rhop)
  teqi0 = getteqi(rho)
  teqim = getteqi(rhom)
  gamp  = 1d0+(alog(teqip)-alog(teqi0))/(alog(rhop)-alog(rho))
  gamm  = 1d0+(alog(teqi0)-alog(teqim))/(alog(rho)-alog(rhom))
  dgdn  = (gamp-gamm)/(rho*epsrho)
  return,dgdn
end

