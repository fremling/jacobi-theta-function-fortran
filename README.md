# Jacobi Theta functions for fortran
Fortran code for jacobi theta funcitons.

The code implements 

$$\theta_3(z,\tau) = \sum_{n=-\infty}^{\infty} e^{\imath \pi \tau n^2} e^{\imath 2 \pi n z}$$


For and extensive list of properties, see the NIST [DLMF website](https://dlmf.nist.gov/20).
Keep howeer in mind that this uses the convention that ads an extra $\pi$ to z compaared no NIST.
This means that here 
$\theta_3(z+1,\tau) = \theta_3(z,\tau)$
wheas NIST has
$\theta_3(z+\pi,\tau) = \theta_3(z,\tau).$


The main file jacobitheta.f90 is stand alone and can be included using `use jacobitheta`.

Most of the rest of the files exists to test the implementation. Run the test with `make test`.
	
All four basic jacobi theta functions are implemented

    theta1(z,tau)
    theta2(z,tau)
    theta3(z,tau)
    theta4(z,tau)
    
as well as the gernalized jacobi theta function `thetagen(a,b,z,tau)` as 

$$\theta_{a,b}(z,\tau) = \sum_{n=-\infty}^{\infty} e^{\imath \pi \tau (n+a)^2} e^{\imath 2 \pi (n+a) (z+b)}$$

where

    theta1(z,tau)=thetagen(.5,.5,z,tau)
    theta2(z,tau)=thetagen(.5, 0,z,tau)
    theta3(z,tau)=thetagen( 0, 0,z,tau)
    theta4(z,tau)=thetagen( 0,.5,z,tau)

The code automatically implements modular transformations that ensures `im(tau) > 1` and `|re(tau)| < 1/2`.

The code also has the option to compute the logarithm of the theta functions using

    logtheta1(z,tau)
    logtheta2(z,tau)
    logtheta3(z,tau)
    logtheta4(z,tau)
    logthetagen(a,b,z,tau)

and this is actually the form used when performing modular transformations. 

A special tanks to Stéphane Laurent for pointin out a few bugs when working on an [R-implementation](https://github.com/stla/jacobi) of the [theta-funcitons](https://github.com/stla/jacobi).

