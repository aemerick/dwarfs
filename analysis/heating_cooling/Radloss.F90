!
!       The cooling function is from Sarazin (1986) Rev.Mod.Phys.
!       and from Raymond 1976
!       
!       ROUTINE FROM PERES ET AL. 1982, APJ 252, 791

!!****f* source/source_terms/cool/radloss/radloss
!!
!! NAME
!!  
!!  radloss
!!
!!
!! SYNOPSIS
!! 
!!  call radloss(T, radia)
!!  call radloss(real, real)
!! 
!! 
!! DESCRIPTION
!!
!!  12/4/05 MKRJ
!!          modified to follow the Dalgarno & McCray cooling curve 
!!          with ionization fraction of 10^-2,
!!          using nested if statements
!!
!! ARGUMENTS
!!
!!  T: 
!!
!!  radia: 
!!
!!***

        subroutine Radloss(T,radia)
!
        implicit none

        real, intent(IN) :: T
        real, intent(OUT) :: radia

!
!  The factor of 5.647e-101 comes from changing units from cgs to kpc, Myr, & Msun
!
	if (T .ge. 1.70e4) then
	  if (T .ge. 5.62e5) then
	    if (T .ge. 2.75e6) then
	      if (T .ge. 3.16e7) then
			radia = 3.090e-27*sqrt(T) !*5.647d-101
	      else 
			radia = 5.188e-21/T**0.33 !*5.647d-101
	      endif
	    else 
	      if (T .ge. 1.78e6) then
			radia = 3.890e-4/T**2.95 !*5.647d-101
	      else 
			radia = 1.3e-22*(T/1.5e5)**0.01 !*5.647d-101
	      endif
	    endif
	  else 
	    if (T .ge. 7.94e4) then
	      if (T .ge. 2.51e5) then
			radia = 3.98e-11/T**2 !*5.647d-101
	      else 
			radia = 6.31e-22*(T/1.e6)**0.01 !*5.647d-101
	      endif
	    else 
	      if (T .ge. 3.98e4) then
			radia = 1.e-31*T**2 !*5.647d-101
	      else 
			radia = 1.479e-21/(T**0.216d0) !*5.647d-101
	      endif
	    endif
	  endif
	else
	  if (T .ge. 1.e3) then
	    if (T .ge. 6.31e3) then
	      if (T .ge. 1.0e4) then
			radia = 7.63d-81*(T**13.8d0) !*5.647d-101
	      else 
			radia = 3.13e-27*(T**0.396d0) !*5.647d-101
	      endif
	    else 
	      if (T .ge. 3.16e3) then
			radia = 2.64e-26*(T**0.152d0) !*5.647d-101
	      else 
			radia = 5.28e-27*(T**0.352d0) !*5.647d-101
	      endif
	    endif
	  else 
	    if (T .ge. 3.98e1) then
	      if (T .ge. 2.00e2) then
			radia = 3.06e-27*(T**0.431d0) !*5.647d-101
	      else 
			radia = 1.52e-28*(T**0.997d0) !*5.647d-101
	      endif
	    else 
	      if (T .ge. 2.51e1) then
			radia = 2.39e-29*(T**1.50d0) !*5.647d-101
	      else 
			radia = 1.095e-32*(T**3.885d0) !*5.647d-101
	      endif
	    endif
	  endif
        endif
! don't cool below 10 K
			if (T .le. 1.e1) radia = 0.e-28 !*5.647d-101  !below T=0, no cooling

        return
        end
