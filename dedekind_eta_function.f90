module dedekind_eta_function
  USE typedef
  
  implicit none

  complex(kind=dpc) :: tau_stored,eta_stored,log_eta_stored
  logical :: eta_set=.FALSE.
  
  private
  public dedekind_eta
  public log_dedekind_eta
  
contains
  
  complex(KIND=dpc) function dedekind_eta(tau) result(eta)
    complex(KIND=dpc), intent(IN)  :: tau    
    eta = do_dedekind_eta(tau,.false.)
  end function dedekind_eta
  
  complex(KIND=C_DOUBLE_COMPLEX) function log_dedekind_eta(tau) result(eta) bind(c)
    use ISO_C_BINDING
    complex(KIND=C_DOUBLE_COMPLEX), intent(IN)  :: tau    
    eta = do_dedekind_eta(tau,.true.)
  end function log_dedekind_eta
  

  complex(KIND=dpc) function do_dedekind_eta(tau,logged) result(eta)
    !!Comptes the eta_function if needed
    !!If already computed lopoks it up
    complex(KIND=dpc), intent(IN)  :: tau
    logical, intent(IN) :: logged
    
    if(eta_set)then
       !write(*,*) 'Eta already set'
       if(tau_stored.eq.tau)then
          !write(*,*) 'with correct tau'
          !!Set eta to the strored one
          if(logged)then
             eta=log_eta_stored
          else
             eta=eta_stored
          end if
       else
          !write(*,*) 'with wrong tau'
          eta=compute_store_eta(tau,logged)
       end if
    else
       !write(*,*) 'Eta not set'
       eta=compute_store_eta(tau,logged)
    end if
    
  end function do_dedekind_eta
  
  complex(KIND=dpc) function compute_store_eta(tau,logged) result(eta)  
    complex(KIND=dpc), intent(IN)  :: tau
    logical, intent(IN) :: logged
    !!Compute eta and store the value
    !write(*,*) 'Compute and store with tau=',tau
    eta=dedekind_eta_compute(tau) 
    eta_set=.TRUE.
    eta_stored=eta
    tau_stored=tau
    log_eta_stored=log(eta)
    if(logged)then
       eta=log(eta)
    end if
    
  end function compute_store_eta
  

complex(KIND=dpc) function dedekind_eta_compute(tau) result(eta)
    !!Computes the dedekind eta function
    complex(KIND=dpc), intent(IN)  :: tau

    if(aimag(tau).ge.1.d0)then
       eta=do_dedekind_eta_compute(tau)
    else
       eta=sqrt(iunit/tau)*do_dedekind_eta_compute(-1.d0/tau)
    end if

  end function dedekind_eta_compute

  complex(KIND=dpc) function do_dedekind_eta_compute(tau) result(eta)
    !!Computes the dedekind eta function
    complex(KIND=dpc), intent(IN)  :: tau
    complex(KIND=dpc) :: q,q24 !The nome
    complex(KIND=dpc) :: increment !!
    integer :: n !!Incerement index
    
    !!Eta start
    !write(*,*) 'tau in:',tau
    q=exp(iunit*2*pi*tau) !!The nome
    q24=exp(iunit*pi*tau/12) !!A 24th of the nome
    eta=1.d0+0*iunit
    !write(*,*) 'Nome;',q
    !write(*,*) 'Nome/24;',q24
    
    n=1
    do 
       !write(*,*) 'n=',n
       increment = (-1)**n * (q**(((n**2)*3+n)/2)+q**(((n**2)*3-n)/2))
       !write(*,*) 'increment:',increment
       !write(*,*) 'quotient:',increment/eta
       eta = eta+increment
       !write(*,*) 'eta_n:',eta
       !!if((abs(increment/eta)).lt.(1d-30))then
       if((eta-increment).eq.eta)then
          !write(*,*) 'DONE'
          exit !! the quotient is small enought
       end if
       n=n+1
    end do
    eta=eta*q24
    !write(*,*) 'eta out:',eta
    
  end function do_dedekind_eta_compute
  
  
end module dedekind_eta_function
