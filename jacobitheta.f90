module jacobitheta
  !!No dependencies -- this is a stand-alone application
  
  implicit none
  !imaginary unit
  COMPLEX(KIND=KIND((1.0D0,1.0D0))), PARAMETER :: iunit = (0.d0, 1.d0)
  !imagianry zero
  COMPLEX(KIND=KIND((1.0D0,1.0D0))), PARAMETER :: czero = (0.d0, 0.d0)
  !pi
  real(KIND=KIND(1.0D0)), parameter :: pi=3.14159265358979323846264338327950d0
  
  private
  public :: theta1,theta2, theta3, theta4, thetagen
  public :: logtheta1,logtheta2, logtheta3, logtheta4, logthetagen
  
contains

  !! All pubblic functions map onto theta3 
  !! Theta 3 is then computed in logarithmic form (to avoid infinites) 
  !! and then transformed back to normal if needed.

  !! All theta functions are (a)periodic with periodicity '1'

  !!------------------------------------------------
  !!         Externally callable functions
  !!------------------------------------------------


  !! ------------------  !!! Theta 1
  function theta1(z,tau,maxiter) result(res)
    complex(KIND=KIND((1.0D0,1.0D0))), intent(IN) :: z,tau
    complex(KIND=KIND((1.0D0,1.0D0))) :: res
    integer, intent(in), optional :: maxiter !!Loop counter!!
    res = dotheta1(z,tau,.FALSE.,maxiter=maxiter)
  end function theta1
  
  !! ------------------  !!! log Theta 1
  function logtheta1(z,tau,maxiter) result(res) bind(C)
    use ISO_C_BINDING
    complex(KIND=C_DOUBLE_COMPLEX), intent(IN) :: z,tau
    complex(KIND=C_DOUBLE_COMPLEX) :: res
    integer, intent(in), optional :: maxiter !!Loop counter!!
    res = dotheta1(z,tau,.TRUE.,maxiter=maxiter)
  end function logtheta1
  
  !! ------------------  !!! Theta 2
  complex(KIND=KIND((1.0D0,1.0D0))) function theta2(z,tau,maxiter) result(res)
    complex(KIND=KIND((1.0D0,1.0D0))), intent(IN) :: z,tau
    integer, intent(in), optional :: maxiter !!Loop counter!!
    res = dotheta2(z,tau,.FALSE.,maxiter=maxiter)
  end function theta2
  
  !! ------------------  !!! log Theta 2
  function logtheta2(z,tau,maxiter) result(res) bind(C)
    use ISO_C_BINDING
    complex(KIND=C_DOUBLE_COMPLEX), intent(IN) :: z,tau
    complex(KIND=C_DOUBLE_COMPLEX) :: res
    integer, intent(in), optional :: maxiter !!Loop counter!!
    res = dotheta2(z,tau,.TRUE.,maxiter=maxiter)
  end function logtheta2

  !! ------------------  !!! Theta 3
  function theta3(z,tau,maxiter) result(res)
    complex(KIND=KIND((1.0D0,1.0D0))), intent(IN) :: z,tau
    complex(KIND=KIND((1.0D0,1.0D0))) :: res
    integer, intent(in), optional :: maxiter !!Loop counter!!
    res = dotheta3(z,tau,.FALSE.,maxiter=maxiter)
  end function theta3
  
  !! ------------------  !!! log Theta 3
  function logtheta3(z,tau,maxiter) result(res) bind(C)
    use ISO_C_BINDING
    complex(KIND=C_DOUBLE_COMPLEX), intent(IN) :: z,tau
    complex(KIND=C_DOUBLE_COMPLEX) :: res
    integer, intent(in), optional :: maxiter !!Loop counter!!
    res = dotheta3(z,tau,.TRUE.,maxiter=maxiter)
  end function logtheta3

  !! ------------------  !!! Theta 4
  function theta4(z,tau,maxiter) result(res)
    complex(KIND=KIND((1.0D0,1.0D0))), intent(IN) :: z,tau
    complex(KIND=KIND((1.0D0,1.0D0))) :: res
    integer, intent(in), optional :: maxiter !!Loop counter!!
    res = dotheta4(z,tau,.FALSE.,maxiter=maxiter)
  end function theta4
  
  !! ------------------  !!! log Theta 4
  function logtheta4(z,tau,maxiter) result(res) bind(C)
    use ISO_C_BINDING
    complex(KIND=C_DOUBLE_COMPLEX), intent(IN) :: z,tau
    complex(KIND=C_DOUBLE_COMPLEX) :: res
    integer, intent(in), optional :: maxiter !!Loop counter!!
    res = dotheta4(z,tau,.TRUE.,maxiter=maxiter)
  end function logtheta4
  
  !! ------------------  !!! Theta gen
  function thetagen(a,b,z,tau,maxiter) result(res)
    complex(KIND=KIND((1.0D0,1.0D0))), intent(IN) :: z,tau
    real(KIND=KIND(1.0D0)), intent(IN) :: a,b
    complex(KIND=KIND((1.0D0,1.0D0))) :: res
    integer, intent(in), optional :: maxiter !!Loop counter!!
    !! thetagen is periodic in a
    res = dothetagen(mod(a,1.0),z+b,tau,.FALSE.,maxiter=maxiter)
  end function thetagen
  
  function logthetagen(a,b,z,tau,maxiter) result(res)
    complex(KIND=KIND((1.0D0,1.0D0))), intent(IN) :: z,tau
    real(KIND=KIND(1.0D0)), intent(IN) :: a,b
    complex(KIND=KIND((1.0D0,1.0D0))) :: res
    integer, intent(in), optional :: maxiter !!Loop counter!!
    res = dothetagen(mod(a,1.0d0),z+b,tau,.TRUE.,maxiter=maxiter)
  end function logthetagen
  
  
  !!--------------------------------------------------
  !! Internal functions
  !!--------------------------------------------------
  
  function M(z,tau) result(res)
    complex(KIND=KIND((1.0D0,1.0D0))), intent(IN) :: z,tau
    complex(KIND=KIND((1.0D0,1.0D0))) :: res
    res = iunit*pi*z+iunit*pi*tau/4
  end function M
  
  !! NB: Thus theta[a,b](z,tau) = theta[a,0](z+b,tau)
  function dothetagen(a,z,tau,logged,maxiter) result(res)
    !!thetagernal = sum_k exp(i pi tau (k+a)^2) exp(i 2 (k +a) (z+b))
    complex(KIND=KIND((1.0D0,1.0D0))) :: res
    complex(KIND=KIND((1.0D0,1.0D0))), intent(IN) :: z,tau
    real(KIND=KIND(1.0D0)), intent(IN) :: a
    logical, intent(IN) :: logged
    integer, intent(in), optional :: maxiter !!Loop counter!!
    res = iunit*pi*a*(2*z+a*tau) + dologtheta3(z+a*tau,tau,maxiter=maxiter)
    if(.not.logged)then
       res = exp(res)
    end if
  end function dothetagen
  
  function dotheta1(z,tau,logged,maxiter) result(res)
    complex(KIND=KIND((1.0D0,1.0D0))) :: res
    complex(KIND=KIND((1.0D0,1.0D0))), intent(IN) :: z,tau
    logical, intent(IN) :: logged
    integer, intent(in), optional :: maxiter !!Loop counter!!
    res = -iunit*pi*0.5d0 + M(z,tau) + dologtheta3(z+0.5d0+tau/2,tau,maxiter=maxiter)
    if(.not.logged)then
       res = exp(res)
    end if
  end function dotheta1
  
  function dotheta2(z,tau,logged,maxiter) result(res)
    complex(KIND=KIND((1.0D0,1.0D0))) :: res
    complex(KIND=KIND((1.0D0,1.0D0))), intent(IN) :: z,tau
    logical, intent(IN) :: logged
    integer, intent(in), optional :: maxiter !!Loop counter!!
    res = M(z,tau) + dologtheta3(z+tau/2,tau,maxiter=maxiter)
    if(.not.logged)then
       res = exp(res)
    end if
  end function dotheta2
  
  recursive function dotheta3(z,tau,logged,maxiter) result(res)
    complex(KIND=KIND((1.0D0,1.0D0))) :: res
    complex(KIND=KIND((1.0D0,1.0D0))), intent(IN) :: z,tau
    logical, intent(IN) :: logged
    integer, intent(in), optional :: maxiter !!Loop counter!!
    res = dologtheta3(z,tau,maxiter=maxiter)
    if(.not.logged)then
       res = exp(res)
    end if
  end function dotheta3

  recursive function dotheta4(z,tau,logged,maxiter) result(res)
    complex(KIND=KIND((1.0D0,1.0D0))) :: res
    complex(KIND=KIND((1.0D0,1.0D0))), intent(IN) :: z,tau
    logical, intent(IN) :: logged
    integer, intent(in), optional :: maxiter !!Loop counter!!
    res = dologtheta3(z+0.5d0,tau,maxiter=maxiter)
    if(.not.logged)then
       res = exp(res)
    end if
  end function dotheta4
  
  !!Here we do all the modular transformations, everything is in log scale
  recursive function dologtheta4(z,tau,pass_in,maxiter) result(res)
    complex(KIND=KIND((1.0D0,1.0D0))) :: res
    complex(KIND=KIND((1.0D0,1.0D0))), intent(IN) :: z,tau
    integer, intent(in), optional :: pass_in,maxiter !!Loop counter!!
    integer :: passes
    if(.not.present(pass_in))then
       passes=1
    else
       passes=pass_in+1
    end if
    res = dologtheta3(z+0.5d0,tau,passes,maxiter=maxiter)
  end function dologtheta4
  
  recursive function dologtheta3(z,tau,pass_in,maxiter) result(res)
    complex(KIND=KIND((1.0D0,1.0D0))) :: res
    complex(KIND=KIND((1.0D0,1.0D0))), intent(IN) :: z,tau
    complex(KIND=KIND((1.0D0,1.0D0))) :: tau2,tauprime
    integer, intent(in), optional :: pass_in, maxiter !!Loop counter!!
    integer :: passes,local_maxiter
    if(.not.present(maxiter))then
       local_maxiter=1000       
    else
       local_maxiter=maxiter
    end if
    
    if(.not.present(pass_in))then
       passes=1
    else
       passes=pass_in+1
    end if
    if(passes.gt.local_maxiter)then
       !! STDERR = 0 !!
       write(*,*) 'ERROR: More than ',local_maxiter,' modular transformations have been performed!'
       write(*,*) '       This can happen if tau is faar from the fundamental domain!'
       write(*,*) '       Either increase maxiter or find perform some transformations by hand...'
       write(*,*) '       Input z  =',z
       write(*,*) '       Input tau=',tau
       call exit(-1)
       end if

    !!write(*,*) '    passes  =',passes
    !!write(*,*) '    First Input z  =',z
    !!write(*,*) '    FIrst Input tau=',tau

    
    !!Check that computation has a chance of converging
    !! STDERR = 0 !!
    if(Aimag(tau).le.0.d0)then
       write(0,*) 'ERROR JACOBI THETA:'
       write(0,*) 'Attempting to evaluate the function with Im(tau)<=0'
       write(0,*) 'tau=',tau
       write(0,*) 'This can never converge, quitting...'
       call exit(-1)
    end if
    
    !!Initaite to avoid compiler warning     
    tau2=tau

    !! This should put tau in the principal region (approx)

    ! If  |Re(tau)| > 1 shift to |Re(tau)|<1
    if(REAL(tau).ge.(0.d0)) then
       tau2 = mod(real(tau+1),2.d0)-1 + iunit*aimag(tau)
    elseif(REAL(tau).lt.(0.d0)) then
       tau2 = mod(real(tau-1),2.d0)+1 + iunit*aimag(tau)
    end if
    
    ! if |Re(tau)| > .6 shift to |Re(tau)| < .6, here theta3 -> theta4
    if(REAL(tau2).gt.(6.d-1)) then
       res = dologtheta4(z,tau2-1,passes,maxiter=local_maxiter) !! Shift back
    elseif(REAL(tau2).le.(-(6.d-1))) then
       res = dologtheta4(z,tau2+1,passes,maxiter=local_maxiter) !! shift forward
       ! if |tau| < 1 invert to make |tau| > 1
       !!NB: There is a risk that the finite presition will cause both |tau|<1 and |1/tau|<1.
       !!tau =(-4.8751486325802507d-2, 0.99881093935790721d0) is such a case!!
       !!We are then content with tau2<0.98 since this will converet almost as fast anyway.
    elseif((ABS(tau2).lt.0.98d0).AND.(AIMAG(tau2).lt.0.98d0)) then
       !! If tau is small. Invert
       tauprime=-1.0d0/tau2
       res = iunit*pi*tauprime*z**2 &
            + dologtheta3(z*tauprime,tauprime,passes,maxiter=local_maxiter) &
            - log(sqrt(-iunit*tau2)) !! invert
    else
       res = argtheta3(z,tau2,maxiter=local_maxiter)
    end if
    
  end function dologtheta3
  
  recursive function argtheta3(z,tau,pass_in,maxiter) result(res)
    complex(KIND=KIND((1.0D0,1.0D0))) :: res
    complex(KIND=KIND((1.0D0,1.0D0))), intent(IN) :: z,tau
    complex(KIND=KIND((1.0D0,1.0D0))) :: zuse,zmin
    integer, intent(in), optional :: pass_in, maxiter !!Loop counter!!
    integer :: passes,local_maxiter

    integer :: quotient

    if(.not.present(maxiter))then
       local_maxiter=10
    else
       local_maxiter=maxiter
    end if
    
    if(.not.present(pass_in))then
       passes=1
    else
       passes=pass_in+1
    end if

    !write(*,*) '    passes  =',passes
    !write(*,*) '       Input z  =',z
    !write(*,*) '       Input tau=',tau
        
    if(passes.gt.local_maxiter)then
       write(*,*) 'ERROR: More than ',local_maxiter,' shifts of z have been performed!'
       write(*,*) '       This should not happen! Please file a bug report!'
       write(*,*) '       Input z  =',z
       write(*,*) '       Input tau=',tau
       call exit(-1)
    end if

    
    !!Reduce  -0.5 < Re(z) < 0.5
    zuse = mod(real(z),1.d0) + iunit*aimag(z)
    !write(*,*) '----         ',zuse


    !Swith the sign of z if the iunit part is negative
    if(AIMAG(zuse).lt.(-AIMAG(tau)/2)) then    
       res = argtheta3(-zuse,tau,passes,maxiter=local_maxiter)
       !!If the argument is to large, shift it away
    elseif(AIMAG(zuse).ge.(AIMAG(tau)/2)) then
       !!write(*,*) 'The quotient is',AIMAG(zuse)/AIMAG(tau)
       quotient = floor(AIMAG(zuse)/AIMAG(tau)+.5) !!!also automatic flooring
       !!write(*,*) 'The floor is',quotient
       
       zmin = zuse-tau*quotient
       res = -2*pi*quotient*iunit*zmin + argtheta3(zmin,tau,passes,maxiter=local_maxiter) - iunit*pi*tau*quotient*quotient
    else
       res = calctheta3(zuse,tau)
    end if
    
  end function argtheta3
  
  function calctheta3(z,tau) result(res)
    complex(KIND=KIND((1.0D0,1.0D0))) :: res
    complex(KIND=KIND((1.0D0,1.0D0))), intent(IN) :: z,tau
    complex(KIND=KIND((1.0D0,1.0D0))) :: q,qweight
    integer :: n
   
    q = exp(iunit*pi*tau)
    res = (1.d0,0.d0)
    n=0
    do 
       n = n + 1
       !!OBS splitting qweight=2*q**(n**2)*cos(2*n*z) into two expontentials to get around lage numbers
       qweight=exp(iunit*n*pi*(tau*n+2*z))+exp(iunit*n*pi*(tau*n-2*z))
       res = res + qweight
       !! Terminate on convergence (or NaNs or Zeroes)
       if( isnan(abs(res))) then
          write(*,*) 'ERROR: Value is infinite in theta function'
          write(*,*) 'z,tau:',z,tau
          call exit(-1)
       elseif( abs(res).eq.0.d0 ) then
          write(*,*) 'Value is zero'
          exit
       elseif((n.ge.3).and.((res+qweight).eq.res)) then
          !!Converged to within machine precision
          exit
       end if
       
    end do
    
    !!Take log on number (to avoid eventual divergentices in the erlier transformations)
    res=log(res)
    
  end function calctheta3

end module jacobitheta
