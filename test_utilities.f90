module test_utilities

  USE typedef     ! types and definitions

  !!Variables
  IMPLICIT NONE  

  public

  contains

!!!------------------------------------------
!!!      Thes the eqivalence of variables
!!!------------------------------------------

  function test_diff(X1,X2,Acc)
    real(KIND=dp), intent(in) :: X1,X2,Acc
    real(KIND=dp) ::  error
    logical :: test_diff
    
    Error=abs((X1-X2)/X2)
    if(abs(error).gt.Acc)then
       test_diff=.true.
    else
       test_diff=.false.
    end if
  end function test_diff

  function test_diff_cmplx(X1,X2,Acc)
    real(KIND=dp), intent(in) :: Acc
    complex(KIND=dpc), intent(in) :: X1,X2
    real(KIND=dp) ::  error
    logical :: test_diff_cmplx
    
    Error=abs((X1-X2)/X2)
    if(abs(error).gt.Acc)then
       test_diff_cmplx=.true.
    else
       test_diff_cmplx=.false.
    end if
  end function test_diff_cmplx

  logical function test_diff_exp_cmplx(X1,X2,Acc,verbose)
    real(KIND=dp), intent(in) :: Acc
    complex(KIND=dpc), intent(in) :: X1,X2
    real(KIND=dp) ::  error,Eff_Acc
    logical, intent(in), optional :: verbose
    logical :: verb
    
    !!Set verbosity
    if(present(verbose))then
       verb=verbose
    else
       verb=.TRUE.
    end if
    
    Error=abs(log(exp(X1-X2)))
    !!Weigh in the fact that if the sace is differnt differnt number of significant numbers are neede
    Eff_Acc = abs(Acc*ceiling(abs(real(X1,dp))+1))
    if(abs(error).gt.Eff_Acc)then
       test_diff_exp_cmplx=.true.
       if(verb)then
          write(*,*) ' --- TDEC (verbose check) --- '
          write(*,*) '            Error:',abs(error),'>',Acc
          write(*,*) ' Effective  Error:',abs(error),'>',Eff_Acc
          write(*,*) 'Input1:',X1
          write(*,*) 'Input2:',X2
          write(*,*) 'Diff:',log(exp(X1-X2))
       end if
    else
       test_diff_exp_cmplx=.false.

    end if
  end function test_diff_exp_cmplx

  
!!!-------------------------------------
!!!    Complex modifyers
!!!------------------------------------- 
  
  complex(KIND=dpc) function mod_2pi(logz) result(mod_logz)
    !! This function returns the argument of logz to be 
    !! within the range 0 <= Im(logz) < 2*pi.
    complex(KIND=dpc), INTENT(IN) :: logz
    real(KIND=dp) :: im_piece
    im_piece = aimag(logz)/(2*pi)
    im_piece = im_piece - nint(im_piece)
    mod_logz = real(logz,dp) + iunit*2*pi*im_piece
  end function mod_2pi


end module test_utilities
