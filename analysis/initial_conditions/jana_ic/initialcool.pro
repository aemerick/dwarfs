pro initialcool

;rho1 = dblarr(n)
;r = (findgen(n)/n)*2000

G =1d0 ;code ;6.673d-8 ;cgs
rho0 =5d2 ;code ; 1d1 ;cm^-3 ;TOTAL central dens 
;b =1.498d-2 ;code ;1.33d2;*3.0857d18 ;pc->cm
gamm = 5d0/3d0 ; gas has constant gamma
k = 1d0 ;1.38e-16 ;cgs
mu = 1d0

T1 = 7.19494d-1 ;8.102d3  ;found with plotcool.pro ;1 ;code ;1d4 ;K
T2 = 8.13208d1 ;7.52d5;1d2 ;code ;1d6 ;K
cs1 =sqrt(gamm*T1);code, calculate then convert ;sqrt(gamm*k*T1/mu/mh)
cs2 =sqrt(gamm*T2);code, calculate then convert ;sqrt(gamm*k*T2/mu/mh) 
rho2rm = 4.49101d-2;1.512d-2 ;code this is the correct value:1.512d-2 ;code ;1.512d-4 ;cm^-3
rho1rm = 5.13970d0 ;code ;1.407d-2 ;cm^-3
rho2rl =4.49100d-2 ;1d-2  ;this is the background density 1d-2 ;code  ;1d-4 ;in cm^-3
b = 5d2*1.126d-4
cphi = 4d0*!pi*G*rho0*b^3
rl = 3d2*1.126d-4 ;pc ->code
rout = 6.071d3*1.126d-4
cp2 = rho2rl/exp((cphi/(cs2^2*rl))*alog(1d0+(rl/b)))

;plot of rho2 equation to get a better estimate of where the root is

n5=1024*3/4
dr  = rout/double(n5)
rad =(findgen(n5)/n5)*rout + 0.5*dr ; in code units
plotrho = dblarr(n5)
plotrho2= dblarr(n5)
dphidr  = dblarr(n5)
rho2 = 0d0
for i=0,n5-1 do begin
  plotrho2[i] = cp2*exp((cphi/(cs2^2*rad[i]))*alog(1d0+rad[i]/b))
  dphidr[i] = -cphi*(rad[i]/((1d0+rad[i]/b)*b)-alog(1d0+rad[i]/b))/rad[i]^2
  print,format='("i,rad,rl,rho2rl,cp2,rho2 = ",i4,e13.5,e13.5,e13.5,e13.5,e13.5,e13.5)',$
        i,rad[i],rl,rho2rl,cp2,plotrho2[i],dphidr[i]
endfor

for o=0,n5-1 do begin
plotrho[o] =(cp2*exp((cphi/(cs2^2*rad[o]))*alog(1d0+(rad[o]/b))))-rho2rm
;print,cp2,cphi,cs2,rad[o],b,plotrho[o]
endfor
plot,rad/1.126d-4,plotrho/1d2,title='Plot of rho2(r) - rho(rl) to help find Roots & guess rm for Newton Method'

; this is a bisearch for the matching density. The background density has to be
; close to the first-peak density in the cooling curve.

eps = 1d-7
rhi = rout
rlo = dr
rhomid = 1d-32
nmax   = 1000l
ic     = 0l
while((abs(rhomid-rho2rm)/rho2rm gt eps) and (ic le nmax)) do begin
  rmid   = 0.5d0*(rhi+rlo)
  rhomid = cp2*exp((cphi/(cs2^2*rmid))*alog(1d0+rmid/b))
  if (rhomid gt rho2rm) then begin
    rlo = rmid
  endif else begin
    rhi = rmid
  endelse
  ic = ic+1l
  print,'still counting: ',ic,rmid,rhomid,rho2rm
endwhile
if (ic gt nmax) then begin
  print,'no convergence in bisearch. probably rho2rl not close enough to rhopeak'
  stop
endif
print,'hopefully this worked: ',rmid,rhomid

if (1 eq 2) then begin
    ;newton method to find rm
    n = 50
    r = dblarr(n)
    newt = dblarr(n)
    fnewt = dblarr(n)

    ;nn = 100
     ;for j=0,nn-1 do begin
       ;100 initial guesses from 0 to 1 pc
       ;r[0] = 1.126d-4*(j/nn)+1.126d-6 ;code
       ;2.253d-2 ;code  ;2d2;*3.0857d18 ;pc->cm
       r[0]=1d0*1.126d-4 ;pc->code
         for i=1,n-1 do begin
           newt[i] = cp2*exp((cphi/(cs2^2*r[i-1]))*alog(1+(r[i-1]/b)))-rho2rm
           fnewt[i] = cp2*(1+(r[i-1]/b)^(cphi/(cs2^2*r[i-1])))*((cphi/(b*cs2^2*r[i-1]*(1+(r[i-1]/b))))-((cphi*alog(1+(r[i-1]/b)))/(cs2^2*r[i-1]^2)))
           d = newt[i]/fnewt[i]
           r[i] = r[i-1] - d
           print,abs(d)
             if (abs(d) lt 1d-5) then begin
               rm = r[i]
             endif else begin
               rm = 3d2*1.126d-4
             endelse
          endfor
             if (rm eq 3d2*1.126d-4) then begin
               print,'BEWARE: DID NOT CONVERGE, USING A FAKE RM!'
             endif
      ;endfor
    ;endfor
;endfor
endif

rm  = rmid

cp1 = rho1rm/exp((cphi/(cs1^2*rm))*alog(1+(rm/b)))
print,'cs1:',cs1,'cs2:',cs2,'cp1:',cp1,'cp2:',cp2,'rm:',rm

rmarr = dblarr(n5)
plotrho1 = dblarr(n5)
;plotrho2 = dblarr(n5)
rhoarr = dblarr(n5)
rhofinal = dblarr(n5)
temp = dblarr(n5)
press = dblarr(n5)
press1= dblarr(n5)
press2= dblarr(n5)
rhofinal2 = dblarr(n5)
tanh1 = dblarr(n5)
tanh2 = dblarr(n5)
drho1dr= dblarr(n5)
drho2dr= dblarr(n5)
drhodr = dblarr(n5)
accgrav= dblarr(n5)
accprss= dblarr(n5)

;set parameters for roll down
sig0 = 5d-4

for p=0,n5-1 do begin
plotrho1[p]=cp1*exp((cphi/(cs1^2*rad[p]))*alog(1d0+(rad[p]/b)))
; do derivatives
drho1dr[p] = getdrhodr(rad[p],cp1,cs1,cphi,b)
drho2dr[p] = getdrhodr(rad[p],cp2,cs2,cphi,b)
drhodr[p]  = getdrhofdr(rad[p],rm,cp1,cs1,cp2,cs2,sig0,cphi,b)
;  plotrho2[p]=cp2*exp((cphi/(cs2^2*rad[p]))*alog(1d0+(rad[p]/b)))
tanh2[p] =5d-1*(tanh((rad[p]-rm)/sig0)+1d0)
tanh1[p] =1d0-tanh2[p]
;continuous version
rhofinal[p] = (tanh1[p]*plotrho1[p])+(tanh2[p]*plotrho2[p])

if (rad[p] lt rm) then begin
  temp[p]    = T1
  press[p]   = plotrho1[p]*T1
  accprss[p] = -cs1^2*drho1dr[p]
  accgrav[p] = -rhofinal[p]*dphidr[p]
endif else begin
  temp[p]    = T2
  press[p]   = plotrho2[p]*T2
  accprss[p] = -cs2^2*drho2dr[p]
  accgrav[p] = -rhofinal[p]*dphidr[p]
endelse
;temp[p] = ((T2-T1)*tanh2[p])+T1
; total acceleration (should be zero)

acctotal   = accgrav+accprss

press1[p]= plotrho1[p]*T1;cs1^2
press2[p]= plotrho2[p]*T2;cs2^2
;press[p] = rhofinal[p]*temp[p]
;discontinuous transition
if rad[p] lt rm then begin
rhofinal2[p] = plotrho1[p]
endif else begin
rhofinal2[p] = plotrho2[p]
endelse
endfor

;plot,rad,press
;plot,rad,tanh1
;oplot,rad,tanh2;,xrange=[0d0,4d0],yrange=[0d0,1.5d0]
;plot,rad,rhofinal;,xrange=[0d0,1d-1],yrange=[1.4d0,1.5d0]
;plot,rad,rhofinal2,xrange=[0,5d-2],yrange=[-1,5d0]
;plot,alog(rad),alog(rhofinal),psym=4

;proteusread,filename='../cooltest/dens0000',dat=densdat,/single
;proteusread,filename='../cooltest/prss0000',dat=prssdat,/single

;for i=256,1023 do begin
;  print,i,rad[i-256],rhofinal[i-256],densdat[256,i]
;endfor

window,0,xsize=600,ysize=400
minr = min([min(plotrho1),min(plotrho2)])
maxr = max([max(plotrho1),max(plotrho2)])
plot,rad,alog10(rhofinal),yrange=alog10([minr,maxr]),line=0,title='Density'
;oplot,rad,alog10(plotrho1)
;oplot,rad,alog10(rhofinal),line=3
oplot,[rm,rm],alog10([0.1*minr,10.0*maxr]),line=1
;oplot,rad,alog10(densdat[256,256:1023]),line=2
;oplot,
;oplot,rad,tanh2
;oplot,rad,rhofinal
;plot,alog(rad/1.126d-4),alog(rhofinal/1d2),title='Initial Density Profile for dSph Stripping with Cooling',xtitle='ln(Radius (pc))',ytitle='ln(Density (cm^-3))';,xrange=[0,600],yrange=[-5d-3,3d0]
window,1,xsize=600,ysize=400
minp = min([min(press1),min(press2),min(press)])
maxp = max([max(press1),max(press2),max(press)])
plot,rad,alog10(press),yrange=alog10([minp,maxp]),line=0,title='Pressure'
;oplot,rad,alog10(press1)
;oplot,rad,alog10(press2),line=2
oplot,[rm,rm],alog10([0.1*minp,10.0*maxp]),line=1
;oplot,rad,alog10(prssdat[256,256:1023]),line=2

;plot the accelerations
window,2,xsize=600,ysize=400
mina = min([min(accgrav),min(accprss),min(acctotal)])
maxa = max([max(accgrav),max(accprss),max(acctotal)])
plot,rad,accgrav,yrange=[mina,maxa],line=0,title='Acceleration'
oplot,rad,accprss,line=2
oplot,rad,acctotal,line=1
end

; derivatives of rho
function getdrhodr,r,cp,cs,cphi,b
  eps = 1d-3
  dr  = r*eps
  rp  = r+0.5d0*dr
  rm  = r-0.5d0*dr 
  rhop= cp*exp((cphi/(cs^2*rp))*alog(1d0+rp/b))
  rhom= cp*exp((cphi/(cs^2*rm))*alog(1d0+rm/b))
  return,(rhop-rhom)/dr
end

function getdrhofdr,r,rmid,cp1,cs1,cp2,cs2,sig0,cphi,b
  eps    = 1d-3
  dr     = r*eps
  rp     = r+0.5d0*dr
  rm     = r-0.5d0*dr
  tanh2p = 5d-1*(tanh((rp-rmid)/sig0)+1d0)
  tanh2m = 5d-1*(tanh((rm-rmid)/sig0)+1d0)
  tanh1p = 1d0-tanh2p
  tanh1m = 1d0-tanh2m
  rhop   = (tanh1p*cp1*exp((cphi/(cs1^2*rp))*alog(1d0+rp/b)))+(tanh2p*cp2*exp((cphi/(cs2^2*rp))*alog(1d0+rp/b)))
  rhom   = (tanh1m*cp1*exp((cphi/(cs1^2*rm))*alog(1d0+rm/b)))+(tanh2m*cp2*exp((cphi/(cs2^2*rm))*alog(1d0+rm/b)))
  return,(rhop-rhom)/dr
end

function getdlnrhodr,r,cp,cphi,cs,b
  eps = 1d-3
  dr  = r*eps
  rp  = r+0.5d0*dr
  rm  = r-0.5d0*dr
  rhop= alog(cp*exp((cphi/(cs^2*rp))*alog(1d0+rp/b)))
  rhom= alog(cp*exp((cphi/(cs^2*rm))*alog(1d0+rm/b)))
  return,(rhop-rhom)/dr
end
