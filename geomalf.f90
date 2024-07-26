! 
! Copyright (c) 2024, Georgios Panou
! 
! This program is free software: you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by the
! Free Software Foundation, either version 3 of the License, or (at your
! option) any later version.
! 
! This program is distributed in the hope that it will be useful, but WITHOUT 
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License 
! along with this program. If not, see <https://www.gnu.org/licenses/>.
! 
! Authors: Iossifidis C., Koci J. and Panou G. <geopanou@survey.ntua.gr>
! 

program geomalf

implicit none
integer*8 NX, NMAX, LFnum, n, p, q, m, t, n2, np1, nm1, nm2, nm3
real*8 theta, PIHALF, raddeg, cosi_1, cosi_2, coeff, coefn, z, cosix, cosx, sini_1, sini_2, sinix
real*8 A, B, coss, Pi, pnn, pn0, C, Pi1, sins, Tni, DSQRT2, DSQRT3
real*8 anm, bnn, bnm, cnn, cnm, bnn_1, cnn_1
real*8 Pn_2_m, Pn_2_m_2, Pn_m_2, Pnm

parameter (NMAX=1000)

real*8 root(0:2*NMAX+4) ! square root lookup table
real*8 rooti(0:2*NMAX+4) ! square root lookup table

real*8 cosn(0:NMAX) ! multiple angle cosines
real*8 sinn(0:NMAX) ! multiple angle sines

real*8 Pn_0(0:NMAX) ! Legendre Polynomials (m=0)
real*8 Pn_1(0:NMAX) ! Legendre Functions (m=1)
real*8 Pn_even(0:NMAX+2) ! Legendre Functions (n is even)
real*8 Pn_odd(0:NMAX+2) ! Legendre Functions (n is odd)
real*8 Pn(0:NMAX+2) ! Vector of all the Legendre Functions with fixed degree

! Initialize
NX=6 ! The maximum degree, less than NMAX
theta=30.0d0 ! The co-latitude in degrees
PIHALF=atan(1.d0)*2.d0
raddeg=PIHALF/90.d0
theta=raddeg*theta
LFnum=(NX+1)*(NX+2)/2

! Computing the square roots and storing them in the lookup tables
do t=0,2*NMAX+4
    root(t)=dsqrt(1.d0*t)    
enddo
p=int(dsqrt(2.d0*NMAX+4))+1
do t=0,p
    q=t*t
    root(q)=1.d0*t
enddo
DSQRT2=root(2)
DSQRT3=root(3)
coeff=0.0d0
do t=0,2*NMAX+3
    z=root(t+1)
    rooti(t)=coeff*z 
    coeff=z
enddo

! Computes the multiple angle cosines by using the Chebyshev's method
! https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Chebyshev_method

cosi_1=dcos(theta)
coeff=2.d0*cosi_1  
cosi_2=1.d0
cosn(0)=1.d0
cosn(1)=cosi_1
do n=2,NX
    cosix=coeff*cosi_1-cosi_2;
    cosn(n)=cosix
    cosi_2=cosi_1
    cosi_1=cosix
enddo

! Computes the multiple angle sines by using the Chebyshev's method
! https://en.wikipedia.org/wiki/List_of_trigonometric_identities#Chebyshev_method

cosx=2.d0*dcos(theta)
sini_2=0.d0
sini_1=dsin(theta)
sinn(0)=0.d0
sinn(1)=sini_1
do n=2,NX
    sinix=cosx*sini_1-sini_2
    sinn(n)=sinix
    sini_2=sini_1
    sini_1=sinix
enddo

! Computes the normalized Legendre Polynomials (m = 0) and Functions (m = 1)
Pn_0(0)=1.d0
Pn_0(1)=DSQRT3*cosn(1)
Pn_1(0)=0.d0
Pn_1(1)=DSQRT3*sinn(1)

! Computes the odd Legendre Polynomials (m = 0) and Functions (m = 1)
pnn=1.d0
coss=cosn(1)
sins=sinn(1)
do n=3,NX,2
    p=n-1
    A=1.d0/p
    B=1.d0/n
    coeff=1.d0-0.75d0*B
    pnn=pnn*(1.d0-A*coeff)
    Pi=coss
    Pi1=sins
    q=n+2
    t=3
    do while (p.gt.0)
        ! Computes the Tni ratio
        B=1.d0-1.d0/p
        C=1.d0+1.d0/q
        Tni=B*C
        Pi=Pi*Tni
        Pi1=Pi1*Tni
        Pi=Pi+cosn(t)
        Pi1=Pi1+t*sinn(t)
        t=t+2
        p=p-2
        q=q+2
    enddo
    coeff=root(2*n+1)
    Pn_0(n)=Pi*pnn*coeff 
    Pn_1(n)=Pi1*pnn*DSQRT2*coeff/rooti(n)
enddo

! Computes the even Legendre Polynomials (m = 0) and Functions (m = 1)
pn0=1.d0
pnn=2.d0
coss=cosn(2)
sins=sinn(2)
do n=2,NX,2
    p=n-1
    A=1.d0/p
    B=1.d0/n
    coeff=1.d0-0.75d0*B
    pnn=pnn*(1.d0-A*coeff)
    ! Computes the tn ratio
    A=1.d0-B
    pn0=pn0*A*A
    Pi=coss
    Pi1=2.d0*sins
    t=4
    q=n+3
    p=p-1
    do while (p.gt.0)
        ! Calculates the Tni ratio
        B=1.d0-1.d0/p
        C=1.d0+1.d0/q
        Tni=B*C
        Pi=Pi*Tni
        Pi1=Pi1*Tni
        Pi=Pi+cosn(t)
        Pi1=Pi1+t*sinn(t)
        t=t+2
        p=p-2
        q=q+2
    enddo
    coeff=root(2*n+1)
    Pn_0(n)=coeff*(Pi*pnn+pn0)
    Pn_1(n)=Pi1*pnn*DSQRT2*coeff/rooti(n)
enddo

write(*,"(a33,i20)") "N=",NX
write(*,"(a33,1p1e20.8)") "theta(deg)=",theta/raddeg
write(*,"(a33,i20)") "The number of Legendre Functions=",LFnum
write(*,"(2a10,4a27)") "n","m","Pnm"

! Computation of the Legendre Functions for m > 1
n=0
write(*,"(2i10,1p1e27.16)")n,0,Pn_0(0)

n=1
Pn_odd(0)=Pn_0(n)
Pn_odd(1)=Pn_1(n)
write(*,"(2i10,1p1e27.16)")n,0,Pn_odd(0)
write(*,"(2i10,1p1e27.16)")n,1,Pn_odd(1)

n=2
Pn_even(0)=Pn_0(n)
Pn_even(1)=Pn_1(n)
Pn_even(2)=root(3)*(root(5)-Pn_0(2))/3.d0
write(*,"(2i10,1p1e27.16)")n,0,Pn_even(0)
write(*,"(2i10,1p1e27.16)")n,1,Pn_even(1)
write(*,"(2i10,1p1e27.16)")n,2,Pn_even(2)

n=3
Pn_odd(0)=Pn_0(n)
Pn_odd(1)=Pn_1(n)
Pn_odd(2)=root(15)*(root(21)*Pn_0(1)/3.d0-Pn_0(3))/5.d0
Pn_odd(3)=root(15)*(root(14)*Pn_1(1)-Pn_1(3))/15.d0
write(*,"(2i10,1p1e27.16)")n,0,Pn_odd(0)
write(*,"(2i10,1p1e27.16)")n,1,Pn_odd(1)
write(*,"(2i10,1p1e27.16)")n,2,Pn_odd(2)
write(*,"(2i10,1p1e27.16)")n,3,Pn_odd(3)

do n=4,NX

    n2=2*n
    np1=n+1
    nm1=n-1
    nm2=n-2
    nm3=n-3
    coefn=root(n2+1)/root(n2-3)
    
    if (mod(n,2).ne.0) then !========== n is odd ==========!
	
	!*=*=*=*=*=*= m is odd *=*=*=*=*=*=!
        Pn(1)=Pn_1(n)
        do m=3,nm2,2
            z=rooti(nm1+m)
            anm=coefn*rooti(nm1-m)/z
            bnm=coefn*rooti(nm3+m)/z
            cnm=rooti(np1-m)/z
            Pn_m_2=Pn(m-2)
            Pn_2_m_2=Pn_odd(m-2)
            Pn_2_m=Pn_odd(m)
            Pnm=anm*Pn_2_m+bnm*Pn_2_m_2-cnm*Pn_m_2
            Pn(m)=Pnm
        enddo
		
        Pn_2_m_2=Pn_odd(nm2)
        Pn_m_2=Pn(nm2)
        cnn=1.d0/root(n)/root(n2-1) ! cnn coefficient
        bnn=root(n2+1)*root(nm1)*cnn ! bnn coefficient
        Pnm=bnn*Pn_2_m_2-cnn*Pn_m_2 ! Computes the Pnn Legendre Function
        Pn(n)=Pnm
		
	!*=*=*=*=*=*= m is even *=*=*=*=*=*=!
        Pn(0)=Pn_0(n)
        m=2
        z=rooti(nm1+m)
        anm=coefn*rooti(nm1-m)/z
        bnm=DSQRT2*coefn*rooti(nm3+m)/z
        cnm=DSQRT2*rooti(np1-m)/z
        Pn_m_2=Pn(m-2)
        Pn_2_m_2=Pn_odd(m-2)
        Pn_2_m=Pn_odd(m)
        Pnm=anm*Pn_2_m+bnm*Pn_2_m_2-cnm*Pn_m_2
        Pn(m)=Pnm
        do m=4,nm2,2
            z=rooti(nm1+m)
            anm=coefn*rooti(nm1-m)/z
            bnm=coefn*rooti(nm3+m)/z
            cnm=rooti(np1-m)/z
            Pn_m_2=Pn(m-2)
            Pn_2_m_2=Pn_odd(m-2)
            Pn_2_m=Pn_odd(m)
            Pnm=anm*Pn_2_m+bnm*Pn_2_m_2-cnm*Pn_m_2
            Pn(m)=Pnm
        enddo
        Pn_2_m_2=Pn_odd(nm3)
        Pn_m_2=Pn(nm3)
        cnn_1=DSQRT3/root(nm1)/root(n2-1) ! cn,n-1 coefficient 
        bnn_1=cnn_1*root(n2+1)*root(nm2)/DSQRT3 ! bn,n-1 coefficient
        Pnm=bnn_1*Pn_2_m_2-cnn_1*Pn_m_2 ! Computes the Pn,n-1 Legendre function
        Pn(nm1)=Pnm
		
        do m=0,n
            Pnm=Pn(m)
            Pn_odd(m)=Pnm
            write(*,"(2i10,1p4e27.16)")n,m,Pnm !Printing the odd Legendre Functions
        enddo
    
    else !========== n is even ==========!
	
	!*=*=*=*=*=*= m is even *=*=*=*=*=*=!
        Pn(0)=Pn_0(n)
        m=2
        z=rooti(nm1+m)
        anm=coefn*rooti(nm1-m)/z
        bnm=DSQRT2*coefn*rooti(nm3+m)/z
        cnm=DSQRT2*rooti(np1-m)/z
        Pn_m_2=Pn(m-2)
        Pn_2_m_2=Pn_even(m-2)
        Pn_2_m=Pn_even(m)
        Pnm=anm*Pn_2_m+bnm*Pn_2_m_2-cnm*Pn_m_2
        Pn(m)=Pnm
        do m=4,nm2,2
            z=rooti(nm1+m)
            anm=coefn*rooti(nm1-m)/z
            bnm=coefn*rooti(nm3+m)/z
            cnm=rooti(np1-m)/z
            Pn_m_2=Pn(m-2)
            Pn_2_m_2=Pn_even(m-2)
            Pn_2_m=Pn_even(m)
            Pnm=anm*Pn_2_m+bnm*Pn_2_m_2-cnm*Pn_m_2
            Pn(m)=Pnm
        enddo
        Pn_2_m_2=Pn_even(nm2)
        Pn_m_2=Pn(nm2)
        cnn=1.d0/root(n)/root(n2-1) ! cnn coefficient
        bnn=root(n2+1)*root(nm1)*cnn ! bnn coefficient
        Pnm=bnn*Pn_2_m_2-cnn*Pn_m_2 ! Computes the Pnn Legendre function
        Pn(n)=Pnm
		
	!*=*=*=*=*=*= m is odd *=*=*=*=*=*=!
        Pn(1)=Pn_1(n)
        do m=3,nm2,2
            z=rooti(nm1+m)
            anm=coefn*rooti(nm1-m)/z
            bnm=coefn*rooti(nm3+m)/z
            cnm=rooti(np1-m)/z
            Pn_m_2=Pn(m-2)
            Pn_2_m_2=Pn_even(m-2)
            Pn_2_m=Pn_even(m)
            Pnm=anm*Pn_2_m+bnm*Pn_2_m_2-cnm*Pn_m_2
            Pn(m)=Pnm
        enddo
		
        Pn_2_m_2=Pn_even(nm3)
        Pn_m_2=Pn(nm3)
        cnn_1=DSQRT3/root(nm1)/root(n2-1) ! cn,n-1 coefficient 
        bnn_1=cnn_1*root(n2+1)*root(nm2)/DSQRT3 ! bn,n-1 coefficient
        Pnm=bnn_1*Pn_2_m_2-cnn_1*Pn_m_2 ! Computes the Pn,n-1 Legendre function
        Pn(nm1)=Pnm
		
        do m=0,n
            Pnm=Pn(m)
            Pn_even(m)=Pnm
            write(*,"(2i10,1p4e27.16)")n,m,Pnm !Printing the even Legendre Functions
        enddo
    endif
enddo

stop
end program
