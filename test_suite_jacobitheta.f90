program test_suite_jacobitheta
  
  USE typedef     ! types and definitions
  use jacobitheta
  use test_utilities
  use dedekind_eta_function
  use misc_random
  
  !!Variables
  IMPLICIT NONE  
  
  
    !!Program Start here!!
    write(*,*) 'Test the jacobi theta function'
    
    !! Testing jacobi theta
    write(*,*) '-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
    call test_jacobi3_special_tricky_values
    write(*,*) '-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
    call test_jacobi3_real
    write(*,*) '-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
    call test_jacobi3_imag
    write(*,*) '-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
    call test_jacobi3_low_tau
    write(*,*) '-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
    call test_jacobi_gen_a
    write(*,*) '-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
    call test_jacobi3_pbc
    write(*,*) '-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
    call test_jacobi1_apbc
    write(*,*) '-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
    call test_jacobi3_z_tau
    write(*,*) '-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
    call test_jacobi3_z_tau_2
    write(*,*) '-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
    call test_jacobi1_near_zero
    write(*,*) '-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
    call test_dedekind_eta_imag_tau
    write(*,*) '-.-.-.-.-.-.-.-.-.-.-.-.-.-.'
    call test_dedekind_eta_gen_tau
    write(*,*) '-.-.-.-.-.-.-.-.-.-.-.-.-.-.'


contains 

  subroutine test_jacobi3_pbc
    integer :: indx
    complex(KIND=dpc) :: tau,z,res1,res2
    real(KIND=dp) :: x,y
    write(*,*) 'test that jacobi 3 is even under z->z+1'

    CALL INIT_RANDOM_SEED()  !!get new seed every time
    !!Set the center of mass parameters
    
    do indx=1,20 !!Test 20 random confugraitions
       call RANDOM_NUMBER(x)
       call RANDOM_NUMBER(y)
       tau=x-.5d0+iunit*(y+.5D0) !!set tau around tau=i
       call RANDOM_NUMBER(x)
       call RANDOM_NUMBER(y)
       z=x+tau*y !!set z
       res1=theta3(z,tau)
       res2=theta3(z+1.d0,tau)
       
       if(test_diff_cmplx(res1,res2,1.d-12))then
          write(*,*) 'For z=',z
          write(*,*) 'and tau=',tau
          write(*,*) 'periodicity is not true!'
          write(*,*) 'res 1:',res1
          write(*,*) 'res 2:',res2
          write(*,*) 'diff:',res1-res2
          call exit(-2)
       end if
    end do
  end subroutine test_jacobi3_pbc

  subroutine test_jacobi1_apbc
    integer :: indx
    complex(KIND=dpc) :: tau,z,res1,res2
    real(KIND=dp) :: x,y
    write(*,*) 'test that jacobi 1 is odd under z->z+1'
    
    CALL INIT_RANDOM_SEED()  !!get new seed every time
    !!Set the center of mass parameters
    
    do indx=1,20 !!Test 20 random confugraitions
       call RANDOM_NUMBER(x)
       call RANDOM_NUMBER(y)
       tau=x-.5d0+iunit*(y+.5D0) !!set tau around tau=i
       call RANDOM_NUMBER(x)
       call RANDOM_NUMBER(y)
       z=x+tau*y !!set z
       res1=theta1(z,tau)
       res2=-theta1(z+1.d0,tau)
       
       if(test_diff_cmplx(res1,res2,1.d-09))then
          write(*,*) 'For z=',z
          write(*,*) 'and tau=',tau
          write(*,*) 'periodicity is not true!'
          write(*,*) 'res 1:',res1
          write(*,*) 'res 2:',res2
          write(*,*) 'diff: ',res1-res2
          write(*,*) 'realtive_diff: ',abs((res1-res2)/res1)
          call exit(-2)
       end if
    end do
  end subroutine test_jacobi1_apbc


  subroutine test_jacobi3_z_tau
    integer :: indx
    complex(KIND=dpc) :: tau,z,res1,res2
    real(KIND=dp) :: x,y
    write(*,*) 'test that jacobi 3 is invariant under z->z+tau, when z=-tau/2'
    
    
    CALL INIT_RANDOM_SEED()  !!get new seed every time
    !!Set the center of mass parameters
    
    do indx=1,20 !!Test 20 random confugraitions
       call RANDOM_NUMBER(x)
       call RANDOM_NUMBER(y)
       tau=x-.5d0+iunit*(y+.5D0) !!set tau around tau=i
       z=-tau/2 !!set z
       res1=theta3(z,tau)
       res2=theta3(z+tau,tau)
       
       if(test_diff_cmplx(res1,res2,1.d-12))then
          write(*,*) 'For z=',z
          write(*,*) 'and tau=',tau
          write(*,*) 'periodicity is not true!'
          write(*,*) 'res 1:',res1
          write(*,*) 'res 2:',res2
          call exit(-2)
       end if
    end do
  end subroutine test_jacobi3_z_tau

  subroutine test_jacobi3_z_tau_2
    integer :: indx
    complex(KIND=dpc) :: tau,z,res1,res2
    real(KIND=dp) :: x,y
    write(*,*) 'test that jacobi 3 is invariant under z->z+tau, for generic z'
    
    
    CALL INIT_RANDOM_SEED()  !!get new seed every time
    !!Set the center of mass parameters
    
    do indx=1,20 !!Test 20 random confugraitions
       call RANDOM_NUMBER(x)
       call RANDOM_NUMBER(y)
       tau=x-.5d0+iunit*(y+.5D0) !!set tau around tau=i
       z=-tau/2 !!set z
       res1=theta3(z,tau)
       res2=exp(-iunit*pi*tau)*exp(-iunit*2*pi*z)*theta3(z+tau,tau)
       
       if(test_diff_cmplx(res1,res2,1.d-12))then
          write(*,*) 'For z=',z
          write(*,*) 'and tau=',tau
          write(*,*) 'periodicity is not true!'
          write(*,*) 'res 1:',res1
          write(*,*) 'res 2:',res2
          call exit(-2)
       end if
    end do
  end subroutine test_jacobi3_z_tau_2

  subroutine test_jacobi3_real()
    integer :: N,Steps=10
    complex(KIND=dpc) :: tau=(0.d0,1.d0),test,z
    real(kind=dp) :: thetalist(11),Error

    thetalist(1)=1.0864348112133080146d0
    thetalist(2)=1.0699237438336250814d0
    thetalist(3)=1.0267020276347579855d0
    thetalist(4)=0.97328668708831650794d0
    thetalist(5)=0.93008056675858800725d0
    thetalist(6)=0.91357913815611682141d0
    thetalist(7)=0.93008056675858800725d0
    thetalist(8)=0.97328668708831650794d0
    thetalist(9)=1.0267020276347579855d0
    thetalist(10)=1.0699237438336250814d0
    thetalist(11)=1.0864348112133080146d0

    write(*,*) 'Testing real argument theta3'
    do N=0,Steps
       z=1.d0*N/Steps
       test=theta3(z,tau)
       Error=(real(test,kind=dp)-thetalist(N+1))/thetalist(N+1)
       if(abs(error).gt.1.d-15)then
          write(*,*) 'tau:',tau
          write(*,*) 'z:',z
          write(*,*) 'theta3-calc:',real(test,kind=dp)
          write(*,*) 'theta3-raw :',thetalist(N+1)
          write(*,*) 'Error:',Error
          call exit(-1)
       end if
    end do
  end subroutine test_jacobi3_real
  

  subroutine test_jacobi3_imag()
    integer :: N,Steps=10
    complex(KIND=dpc) :: tau=(0.d0,1.d0),z
    real(kind=dp) :: thetalist(11),test

    write(*,*) 'Testing imaginary argument of theta3'
    thetalist(1)=1.0864348112133080146d0
    thetalist(2)=1.1040699485314762371d0
    thetalist(3)=1.1641782302616211721d0
    thetalist(4)=1.2913223115220268668d0
    thetalist(5)=1.5375200463815831453d0
    thetalist(6)=2.0037348984882403346d0
    thetalist(7)=2.8820138107244502994d0
    thetalist(8)=4.5371715100832173800d0
    thetalist(9)=7.6673499601233240664d0
    thetalist(10)=13.630057003345058244d0
    thetalist(11)=25.140854031838732728d0

    do N=1,Steps
       z=iunit*N/Steps
       test=real(theta3(z,tau),dp)
       if(test_diff(test,thetalist(N+1),1.d-15))then
          write(*,*) 'tau,z:',tau,z
          write(*,*) 'theta3-calc:',real(test,kind=dp)
          write(*,*) 'theta3-raw :',thetalist(N+1)
          call exit(-2)
       end if
    end do
  end subroutine test_jacobi3_imag



  subroutine test_jacobi3_low_tau()
    integer :: N,Steps=20
    complex(KIND=dpc) :: tau=(0.d0,1.d0),test
    real(kind=dp) :: thetalist(20),Error

    thetalist(1)=1.0864348112133080146d0
    thetalist(2)=1.4194954880837661234d0
    thetalist(3)=1.7323303588980335702d0
    thetalist(4)=2.0000139493694248360d0
    thetalist(5)=2.2360686514584039042d0
    thetalist(6)=2.4494897746873515544d0
    thetalist(7)=2.6457513125537614827d0
    thetalist(8)=2.8284271248149862514d0
    thetalist(9)=3.0000000000031532911d0
    thetalist(10)=3.1622776601685229690d0
    thetalist(11)=3.3166247903554063592d0
    thetalist(12)=3.4641016151377548809d0
    thetalist(13)=3.6055512754639893063d0
    thetalist(14)=3.7416573867739413862d0
    thetalist(15)=3.8729833462074168852d0
    thetalist(16)=4.0000000000000000000d0
    thetalist(17)=4.1231056256176605498d0
    thetalist(18)=4.2426406871192851464d0
    thetalist(19)=4.3588989435406735522d0
    thetalist(20)=4.4721359549995793928d0

    write(*,*) 'Testing small arguments tau arguments of theta3'
    do N=1,Steps
       tau=iunit/N
       test=theta3(czero,tau)
       Error=(real(test,kind=dp)-thetalist(N))/thetalist(N)
       if(abs(error).gt.1.d-15)then
          write(*,*) 'tau:',tau
          write(*,*) 'theta3-calc:',real(test,kind=dp)
          write(*,*) 'theta3-raw :',thetalist(N)
          write(*,*) 'Error:',Error
          call exit(-1)
       end if
    end do
  end subroutine test_jacobi3_low_tau


  subroutine test_jacobi_gen_a()
    integer :: K,N,  Steps=10,Number
    complex(KIND=dpc) :: tau,test
    real(kind=dp) :: thetalist(65),Error,a

    write(*,*) 'Testing "imaginary" argument of generalized theta'
    thetalist(1)=1.0864348112133080146d0 !!tau=i
    thetalist(2)=1.0864348112133080146d0
    thetalist(3)=1.0037348854877390910d0 !!tau=2i
    thetalist(4)=0.41576060259602703231d0
    thetalist(5)=1.0037348854877390910d0
    thetalist(6)=1.0001613990351406940d0 !!tau=3i
    thetalist(7)=0.36608447993144643810d0
    thetalist(8)=0.36608447993144643810d0
    thetalist(9)=1.0001613990351406940d0
    thetalist(10)=1.0000069746847124180d0 !!tau=4i
    thetalist(11)=0.45678956907805841070d0
    thetalist(12)=0.086427836528595596584d0
    thetalist(13)=0.45678956907805841070d0
    thetalist(14)=1.0000069746847124180d0
    thetalist(15)=1.0000003014034550780d0!!tau=5i
    thetalist(16)=0.53353114347282167664d0
    thetalist(17)=0.084503031554652736425d0
    thetalist(18)=0.084503031554652736425d0
    thetalist(19)=0.53353114347282167664d0
    thetalist(20)=1.0000003014034550780d0
    thetalist(21)=1.0000000130248242722d0!!tau=6i
    thetalist(22)=0.59238691304436208204d0
    thetalist(23)=0.12337467676577213081d0
    thetalist(24)=0.017966582042258856541d0
    thetalist(25)=0.12337467676577213081d0
    thetalist(26)=0.59238691304436208204d0
    thetalist(27)=1.0000000130248242722d0
    thetalist(28)=1.0000000005628536915d0!!tau=7i
    thetalist(29)=0.63839453087357880243d0
    thetalist(30)=0.16610833121572574059d0
    thetalist(31)=0.018372793906149352591d0
    thetalist(32)=0.018372793906149352591d0
    thetalist(33)=0.16610833121572574059d0
    thetalist(34)=0.63839453087357880243d0
    thetalist(35)=1.0000000005628536915d0
    thetalist(36)=1.0000000000243231134d0!!tau=8i
    thetalist(37)=0.67523191105318101580d0
    thetalist(38)=0.20788030129801351616d0
    thetalist(39)=0.029233907312429048192d0
    thetalist(40)=0.0037348854634159776289d0
    thetalist(41)=0.029233907312429048192d0
    thetalist(42)=0.20788030129801351616d0
    thetalist(43)=0.67523191105318101580d0
    thetalist(44)=1.0000000000243231134d0
    thetalist(45)=1.0000000000010510970d0!!tau=9i
    thetalist(46)=0.70534668157882875992d0
    thetalist(47)=0.24752015872562547739d0
    thetalist(48)=0.043217405606128458770d0
    thetalist(49)=0.0039157540904684009573d0
    thetalist(50)=0.0039157540904684009573d0
    thetalist(51)=0.043217405606128458770d0
    thetalist(52)=0.24752015872562547739d0
    thetalist(53)=0.70534668157882875992d0
    thetalist(54)=1.0000000000010510970d0
    thetalist(55)=1.0000000000000454220d0!!tau=10i
    thetalist(56)=0.73040269105752847560d0
    thetalist(57)=0.28460954518952392420d0
    thetalist(58)=0.059164717620983799127d0
    thetalist(59)=0.0065736730122758083340d0
    thetalist(60)=0.00077640640785353249447d0
    thetalist(61)=0.0065736730122758083340d0
    thetalist(62)=0.059164717620983799127d0
    thetalist(63)=0.28460954518952392420d0
    thetalist(64)=0.73040269105752847560d0
    thetalist(65)=1.0000000000000454220d0

    Number=1
    do N=1,Steps
       do K=0,N
          a=1.d0*K/N
          tau=iunit*N
          test=thetagen(a,0.d0,iunit*0.d0,tau)
          Error=(real(test,kind=dp)-thetalist(Number))/thetalist(Number)
          if(abs(error).gt.1.d-14)then
             write(*,*) 'tau,a,N,K,num:',tau,a,N,K,Number
             write(*,*) 'theta_gen-calc:',real(test,kind=dp)
             write(*,*) 'theta_gen-raw :',thetalist(Number)
             write(*,*) 'Error:',Error
             call exit(-2)
          end if
          Number=Number+1
       end do
    end do
  end subroutine test_jacobi_gen_a
  

  subroutine test_jacobi1_near_zero()
    !!Test thata theta_1(0|tau)=0
    complex(kind=dpc) :: res,res_log,tau
    integer :: n,n_max=100
    real(kind=dp) :: tau_1,tau_2
    
    write(*,*) 'Test that theta_1(0|tau)=0'
    call init_random_seed
    do n=1,n_max
       call random_number(tau_1)
       call random_number(tau_2)
       tau_1=tau_1-0.5
       tau_2=tau_2+0.5
       tau=tau_1+iunit*tau_2
       
       res = theta1(czero,tau)
       res_log=logtheta1(czero,tau)
       if(test_diff_cmplx(res,czero,1.d-15).or.&
            test_diff_exp_cmplx(exp(res_log),czero,1.d-15))then
          write(*,*) 'tau:',tau
          write(*,*) 'theta(0)     :',res
          write(*,*) 'log(theta(0)):',res_log
          write(*,*) 'ERROR: Not zero'
          call exit(-1)
       end if
       
       
    end do
    
  end subroutine test_jacobi1_near_zero
  
  
  subroutine test_dedekind_eta_imag_tau()
    integer :: N,Steps=10,index
    complex(kind=dpc) :: eta_list(21),tau,eta_facit,eta_compute
    complex(kind=dpc) :: log_eta_facit,log_eta_compute
    write(*,*) 'Test the dedekint eta function with imaginalry tau'

    
    eta_list(1)=(0.49083632513438365970d0,0.d0)
    eta_list(2)=(0.52522802204669445100d0,0.d0)
    eta_list(3)=(0.55841877906043821728d0,0.d0)
    eta_list(4)=(0.59025394191255272032d0,0.d0)
    eta_list(5)=(0.62061754903949295672d0,0.d0)
    eta_list(6)=(0.64942630081274949031d0,0.d0)
    eta_list(7)=(0.67662115500527686852d0,0.d0)
    eta_list(8)=(0.70215608280148734365d0,0.d0)
    eta_list(9)=(0.72598390486112827388d0,0.d0)
    eta_list(10)=(0.74803996263720069675d0,0.d0)
    eta_list(11)=(0.76822542232605665900d0,0.d0)
    eta_list(12)=(0.78639279165469003545d0,0.d0)
    eta_list(13)=(0.80233629864351675501d0,0.d0)
    eta_list(14)=(0.81578898073872370893d0,0.d0)
    eta_list(15)=(0.82642694495296496979d0,0.d0)
    eta_list(16)=(0.83387987650906921032d0,0.d0)
    eta_list(17)=(0.83774606470719159308d0,0.d0)
    eta_list(18)=(0.83761021439735304861d0,0.d0)
    eta_list(19)=(0.83306292622904557611d0,0.d0)
    eta_list(20)=(0.82372150713672994164d0,0.d0)
    eta_list(21)=(0.80925228968134227482d0,0.d0)
    
    write(*,*) 'Testing some arguments to dedekind_eta'
    index=1
    do N=-Steps,Steps
       tau=iunit*exp(-real(N,dp)/Steps)
       eta_facit=eta_list(index)
       eta_compute=dedekind_eta(tau)
       log_eta_facit=log(eta_list(index))
       log_eta_compute=log_dedekind_eta(tau)
       
       
       if(test_diff_cmplx(eta_facit,eta_compute,1.d-13).or.&
            test_diff_cmplx(log_eta_facit,log_eta_compute,1.d-13))then
          write(*,*) 'ERROR:'
          write(*,*) 'Grid point n=',n
          write(*,*) 'tau:',tau
          write(*,*) 'eta_calc: ', eta_compute
          write(*,*) 'eta_facit:',eta_facit
          write(*,*) 'log_eta_calc: ', log_eta_compute
          write(*,*) 'log_eta_facit:',log_eta_facit
          write(*,*) 'e**i*pi*tau/12:',exp(iunit*pi*tau/12)
          write(*,*) 'DIff:    ', eta_compute-eta_facit
          write(*,*) 'DIff log:', log_eta_compute-log_eta_facit
          call exit(-1)
       end if
       index=index+1
    end do
  end subroutine test_dedekind_eta_imag_tau
  

  subroutine test_dedekind_eta_gen_tau()
    integer :: R,N,  Steps=3,index
    complex(kind=dpc) :: eta_list(49),tau,eta_facit,eta_compute
    complex(kind=dpc) :: log_eta_facit,log_eta_compute
    write(*,*) 'Test the dedekint eta function with general tau'

    eta_list(1)=(0.47411148292811931308d0,-0.12703778897291158971d0) 
    eta_list(2)=(0.48337944335425864982d0,-0.08523285428748957200d0) 
    eta_list(3)=(0.48896857419116250815d0,-0.04277919074651812800d0)
    eta_list(4)=(0.49083632513438365970d0,-0.d0)
    eta_list(5)=(0.48896857419116250815d0,+0.04277919074651812800d0) 
    eta_list(6)=(0.48337944335425864982d0,+0.08523285428748957200d0) 
    eta_list(7)=(0.47411148292811931308d0,+0.12703778897291158971d0) 
    eta_list(8)=(0.58008034760221657631d0,-0.15543206068517887846d0) 
    eta_list(9)=(0.59142362248287680822d0,-0.10428649930967291490d0) 
    eta_list(10)=(0.59826268340206099964d0,-0.05233867420560174521d0)
    eta_list(11)=(0.60054336659657605449d0,-0.d0) 
    eta_list(12)=(0.59826268340206099964d0,+0.05233867420560174521d0) 
    eta_list(13)=(0.59142362248287680822d0,+0.10428649930967291490d0) 
    eta_list(14)=(0.58008034760221657631d0,+0.15543206068517887846d0) 
    eta_list(15)=(0.67018985389265717667d0,-0.17957683012606960888d0) 
    eta_list(16)=(0.68343392947144673148d0,-0.12060271635882839391d0) 
    eta_list(17)=(0.69136076292345412727d0,-0.06039243856485591590d0)
    eta_list(18)=(0.69383159208758175087d0,-0.d0)
    eta_list(19)=(0.69136076292345412727d0,+0.06039243856485591590d0) 
    eta_list(20)=(0.68343392947144673148d0,+0.12060271635882839391d0) 
    eta_list(21)=(0.67018985389265717667d0,+0.17957683012606960888d0) 
    eta_list(22)=(0.74204877583656472634d0,-0.19883137022991071905d0) 
    eta_list(23)=(0.75846577837802827279d0,-0.13499956608052770407d0) 
    eta_list(24)=(0.76756214247399513466d0,-0.06590582108774201802d0)
    eta_list(25)=(0.76822542232605665900d0,-0.d0)
    eta_list(26)=(0.76756214247399513466d0,+0.06590582108774201802d0) 
    eta_list(27)=(0.75846577837802827279d0,+0.13499956608052770407d0) 
    eta_list(28)=(0.74204877583656472634d0,+0.19883137022991071905d0) 
    eta_list(29)=(0.79173576249299643907d0,-0.21214495817883756625d0) 
    eta_list(30)=(0.81957140161427703469d0,-0.15250413498425780255d0) 
    eta_list(31)=(0.83111615789797191373d0,-0.06481300496197290143d0)
    eta_list(32)=(0.81966517608781404600d0,-0.d0) 
    eta_list(33)=(0.83111615789797191373d0,+0.06481300496197290143d0) 
    eta_list(34)=(0.81957140161427703469d0,+0.15250413498425780255d0) 
    eta_list(35)=(0.79173576249299643907d0,+0.21214495817883756625d0) 
    eta_list(36)=(0.80956734066191125513d0,-0.21692291514897072897d0) 
    eta_list(37)=(0.87371107547145868148d0,-0.18338229589627030051d0) 
    eta_list(38)=(0.89140292322525982105d0,-0.04899937055862807979d0)
    eta_list(39)=(0.83812578422521199946d0,-0.d0) 
    eta_list(40)=(0.89140292322525982105d0,+0.04899937055862807979d0) 
    eta_list(41)=(0.87371107547145868148d0,+0.18338229589627030051d0) 
    eta_list(42)=(0.80956734066191125513d0,+0.21692291514897072897d0) 
    eta_list(43)=(0.78167768658677098652d0,-0.20944990486235352356d0) 
    eta_list(44)=(0.93090134364302308322d0,-0.23546353688284172806d0) 
    eta_list(45)=(0.96012409732447025159d0,-0.01349468542183536813d0)
    eta_list(46)=(0.80925228968134227482d0,-0.d0)
    eta_list(47)=(0.96012409732447025159d0,+0.01349468542183536813d0) 
    eta_list(48)=(0.93090134364302308322d0,+0.23546353688284172806d0) 
    eta_list(49)=(0.78167768658677098652d0,+0.20944990486235352356d0)

    write(*,*) 'Testing some arguments to dedekind_eta'
    index=1
    do N=-Steps,Steps
       do R=-Steps,Steps
          tau=iunit*exp(-real(N,dp)/Steps)+real(R,dp)/Steps
          eta_facit=eta_list(index)
          eta_compute=dedekind_eta(tau)

          log_eta_facit=log(eta_list(index))
          log_eta_compute=log_dedekind_eta(tau)
          
       
          if(test_diff_cmplx(eta_facit,eta_compute,1.d-13).or.&
               test_diff_cmplx(log_eta_facit,log_eta_compute,1.d-13))then
             write(*,*) 'ERROR:'
             write(*,*) 'Grid point (n,r):', n,r
             write(*,*) 'tau:',tau
             write(*,*) 'eta_calc:  ', eta_compute
             write(*,*) 'eta_facit :',eta_facit
             write(*,*) 'log_eta_calc:  ', log_eta_compute
             write(*,*) 'log_eta_facit :',log_eta_facit
             write(*,*) 'e**i*pi*tau/12:',exp(iunit*pi*tau/12)
             write(*,*) 'DIff:    ', eta_compute-eta_facit
             write(*,*) 'DIff log:', log_eta_compute-log_eta_facit
             call exit(-1)
          end if
          index = index +1
       end do
    end do
  end subroutine test_dedekind_eta_gen_tau




  subroutine test_jacobi3_special_tricky_values()
    integer :: N
    integer, parameter :: Steps=5
    complex(KIND=dpc) :: test,z,tau
    complex(KIND=dpc) :: taulist(Steps),zlist(Steps),thetalist(Steps)
    real(kind=dp) :: Error

    !!Entry No.1
    !!This caused the jacobi-theta-to go into a loop.
    thetalist(1)=(0.82576066577152429d0,0.13279039838532530d0)
    taulist(1)=(-0.50d0,0.525d0)
    zlist(1)=(0.3d0,7.2d-2)
    !!Entry No.2
    !!This caused the jacobi-theta-to go into a loop.
    thetalist(2)=(0.96477107519937078d0,-3.3912040458649421d-2)
    taulist(2)=(-4.8751486325802507d-2, 0.99881093935790721d0)
    zlist(2)=(0.3d0,7.25d-2)
    !!Entry No.3
    !!This caused the jacobi-theta-to go into a loop.
    thetalist(3)=(-17.861938755626756d0     ,  148.48032551299363d0)
    taulist(3)=(0.37202495190115092d0     , 0.92822272928588978d0)
    zlist(3)=(-0.15500984054594991d0     ,  1.2618627643852796d0)


    thetalist(4)=(655869.94733675814d0     ,  592671.58801614877d0)
    taulist(4)=(0.16552260224358278d0     , 1.01d0)
    zlist(4)=(-0.11073905127262d0     ,  -2.1d0)

    
    thetalist(5)=(-46.426517230164478d0,-41.380842988230604d0)
    taulist(5)=(0.7792256d0     , 1.0d-7)
    zlist(5)=(-0.0d0     ,  0.0d0)


    write(*,*) 'Testing special tricky values for theta3'
    do N=1,Steps
       z=zlist(N)
       tau=taulist(N)
       test=theta3(z,tau)
       Error=abs(test-thetalist(N))
       if(abs(error).gt.1.d-15)then
          write(*,*) 'tau:',tau
          write(*,*) 'z:',z
          write(*,*) 'theta3-calc:',test
          write(*,*) 'theta3-raw :',thetalist(N)
          write(*,*) 'Error:',Error
          call exit(-1)
       end if
    end do
  end subroutine test_jacobi3_special_tricky_values
  


end program test_suite_jacobitheta


