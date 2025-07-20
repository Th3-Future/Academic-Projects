!****************Begin Header***************************************************
!This program was written by Godswill Ezeorah, Student Number: 501012886 on July 20, 2021.
!This program solves an unsteady linear/non-linear equation using finite difference method
!and was written as a solution to AE8112 PS5 q1(b,c,d)
!****************************End Header*****************************************
program unsteady_Dfinite
    implicit none
    !Variable declaration
    DOUBLE PRECISION, dimension(:), ALLOCATABLE :: omg_i, V_r, V_z, r_i, z_i
    DOUBLE PRECISION, PARAMETER :: rh=8.3D-4, RR=1, L=10, mu=1.4D-4, CLVZ=5
    DOUBLE PRECISION, PARAMETER :: Pi = 4*atan(1.0)
    DOUBLE PRECISION :: dr, dz, dt, tols, rrr, start, finish
    INTEGER :: nr, nz, jj, i1, j1, l1, nnr, nnz
    
    call cpu_time(start) !gets the start time, for timing purpose
    open(2, file = 'PS5_Q1w.txt', status = 'unknown') 
    open(3, file = 'PS5_Q1vr.txt', status = 'unknown') 
    open(4, file = 'PS5_Q1vz.txt', status = 'unknown')  
    Nr=21   !Number of grid points along radial direction
    Nz=201  !Number of grid points along axial direction
    nnr=nr-2; nnz=nz-2 !Number of grid points in the interior region
    dr=RR/(nr-1)
    dz=L/(nz-1)
    dt=1D-6     !given timestep
    tols=1D-5
    ALLOCATE(omg_i(nr*nz), V_r(nr*nz), V_z(nr*nz), r_i(nr), z_i(nz))

    !Initialization
    omg_i=0; V_r=0; V_z=0; r_i=0; z_i=0 !first initial consideration
    !rrr=0.75; omg_i=2*CLVZ*(rrr/RR**2); V_r=0; V_z=CLVZ*(1-(rrr/1)**2) !second initial consideration

    call unsteady(omg_i, V_r, V_z, r_i, z_i) !solves the unsteady system

    !Printing and writing results
    j1=nr*nz-(nr-1)
    l1=nr*nz
    jj=nz
    do i1 = j1, 1, -nr
        write(2,*) z_i(jj), omg_i(i1:l1)
        write(3,*) z_i(jj), V_r(i1:l1)
        write(4,*) z_i(jj), V_z(i1:l1)
        l1=l1-nr
        jj=jj-1
    end do
    call cpu_time(finish) !gets the end time
    write(2,*) r_i(1:nr); write(3,*) r_i(1:nr); write(4,*) r_i(1:nr)
    print *, "Computation Time = ", (finish-start)/(60*60),"hours"
    close(2);close(3);close(4)
    
    contains
!*************************************************************************************
    subroutine unsteady(oma, Vr, Vz, ri, zi)
        DOUBLE PRECISION, dimension(nr*nz) :: oma, Vz, omo, b, Vr
        DOUBLE PRECISION, dimension(nnr*nnz):: d, e, f, g, h, bi, Vri
        DOUBLE PRECISION, dimension(nr) :: ri
        DOUBLE PRECISION, dimension(nz) :: zi, mdot, som
        DOUBLE PRECISION :: aa,aw,as,ap,ae,an, k1, k2, k3, t, rp, zp
        DOUBLE PRECISION :: Mmdot, Amdot
        integer  :: i, j, k, nj, nk, ll, ij, iii
    !!This subroutine solves unsteady fluid flow through a pipe problem using vorticity method

        open(1, file = 'PS5_Q1.txt', status = 'unknown')
        !Variable initialization
        d=0; e=0; f=0; g=0; h=0
        t=1D-8; rp=0; zp=0
        b=0
        Vri=0
        Mmdot=1
        ij=1
        do i = 1,nz !Loop for creating the axial and radial gridpoint dimensions
            if ( i<=nr ) then
                ri(i)=rp
                zi(i)=zp
            else
                zi(i)=zp
            end if
            rp=rp+dr
            zp=zp+dz
        end do
        do while (Mmdot>tols) !Main loop for time step
            j=nr+1
            k=nr+nr
            ll=1
            omo=oma

            !The below loop solves the vorticity at each gridpoint
            do i = 1,nr*nz 
                aw=mu*dt/(rh*dr**2)-mu*dt/(2*rh*ri(ll)*dr)+Vr(i)*dt/(2*dr)
                as=mu*dt/(rh*dz**2)+Vz(i)*dt/(2*dz)
                ap=2*mu*dt/(rh*dr**2)+2*mu*dt/(rh*dz**2)+mu*dt/(rh*ri(ll)**2)-Vr(i)*dt/ri(ll) -1
                ae=mu*dt/(rh*dr**2)+mu*dt/(2*rh*ri(ll)*dr)-Vr(i)*dt/(2*dr)
                an=mu*dt/(rh*dz**2)-Vz(i)*dt/(2*dz)
                k2=(Vr(i+nr)/dz)-((Vz(i+1)-Vz(i-1))/(2*dr))
                k3=Vz(i-1)/(2*dr)
                !The 9 if-statements below are for the 9 regions of our CV formulation
                if ( i==1 ) then !bottom left region
                    oma(i)=0 !symmetry condition
                    ll=ll+1
                else if ( i < nr ) then        !bottom region
                    oma(i)=k2
                    ll=ll+1
                else if ( i==nr ) then         !bottom right region
                    oma(i)=k3
                    ll=1
                else if ( i==nr*nz-(nr-1) ) then !top left region
                    oma(i)=0 !symmetry condition
                    ll=ll+1
                else if ( i==nr*nz ) then       !top right region
                    oma(i)=k3 !wall condition
                else if ( i==j ) then         !left region
                    oma(i)=0 !symmetry condition
                    j=j+nr
                    ll=ll+1
                else if ( i==k ) then         !right region
                    oma(i)=k3 !wall condition
                    k=k+nr
                    ll=1
                else if ( i < (nr*nz)-nr ) then    !Interior region
                    oma(i)=aw*omo(i-1)+as*omo(i-nr)-ap*omo(i)+ae*omo(i+1)+an*omo(i+nr)
                    ll=ll+1
                else                          !top region
                    oma(i)=oma(i-nr) !outflow condition
                    ll=ll+1
                end if 
            end do
            !print *,'Oma=', oma(1:nr*nz)

            !The below Loop composes the radial velocity at the boundary
            j=nr+1
            k=nr+nr
            ll=1
            do i = 1,nr*nz
                !The 9 if-statements below are for the 9 regions of our formulation
                if ( i==1 ) then !bottom left region
                    Vr(i)=0 !inflow or symmetry condition
                    b(i)=(ri(ll)**2)*(0-0)/(2*dz)
                    ll=ll+1
                else if ( i < nr ) then        !bottom region
                    Vr(i)=0 !inflow condition
                    b(i)=(ri(ll)**2)*(oma(i+nr)+((Vz(i+1)-Vz(i-1))/(2*dr)))/(2*dz)
                    ll=ll+1
                else if ( i==nr ) then         !bottom right region
                    Vr(i)=0 !wall conditions
                    b(i)=(ri(ll)**2)*(oma(i+nr)-(Vz(i-1)/(2*dr)))/(2*dz)
                    ll=1
                else if ( i==nr*nz-(nr-1) ) then !top left region
                    Vr(i)=0 !symmetry condition
                    b(i)=ri(ll)**2*(0-0)/(2*dz)
                    ll=ll+1
                else if ( i==nr*nz ) then       !top right region
                    Vr(i)=0
                    b(i)=(ri(ll)**2)*(oma(i)-oma(i-nr-nr))/(2*dz)
                else if ( i==j ) then         !left region
                    Vr(i)=0 !symmetry conditions
                    b(i)=ri(ll)**2*(0-0)/(2*dz)
                    j=j+nr
                    ll=ll+1
                else if ( i==k ) then         !right region
                    Vr(i)=0 !wall conditions
                    b(i)=(ri(ll)**2)*(oma(i+nr)-oma(i-nr))/(2*dz)
                    k=k+nr
                    ll=1
                else if ( i < (nr*nz)-nr ) then    !Interior region
                    b(i)=(ri(ll)**2)*(oma(i+nr)-oma(i-nr))/(2*dz)
                    ll=ll+1
                else                          !top region
                    b(i)=(ri(ll)**2)*(oma(i)-oma(i-nr-nr))/(2*dz)
                    ll=ll+1
                end if  
            end do
            nj=nnr+1; nk=nnr+nnr
            j=nr+1; k=1
            ll=2
            do i = 1,nnr*nnz !This loop composes the interior radial velocity pentadiagonal vectors
                aw=(ri(ll)**2)/(dr**2) - ri(ll)/(2*dr)
                aa=(ri(ll)**2)/(dz**2)
                ap=(2*ri(ll)**2)/(dr**2) + (2*ri(ll)**2)/(dz**2) + 1
                ae=(ri(ll)**2)/(dr**2) + ri(ll)/(2*dr)
                !The 9 if-statements below are for the 9 regions of our formulation
                if ( i==1 ) then !bottom left region
                    f(i)=-ap
                    g(i+1)=ae
                    h(i+nnr)=aa
                    bi(i)=b(j+k)
                    k=k+1
                    ll=ll+1
                else if ( i < nnr ) then        !bottom region
                    e(i-1)=aw
                    f(i)=-ap
                    g(i+1)=ae
                    h(i+nnr)=aa
                    bi(i)=b(j+k)
                    k=k+1
                    ll=ll+1
                else if ( i==nnr ) then         !bottom right region
                    e(i-1)=aw
                    f(i)=-ap
                    h(i+nnr)=aa
                    bi(i)=b(j+k)
                    j=j+nr
                    k=1
                    ll=2
                else if ( i==nnr*nnz-(nnr-1) ) then !top left region
                    d(i-nnr)=aa
                    f(i)=-ap+aa
                    g(i+1)=ae
                    bi(i)=b(j+k)
                    k=k+1
                    ll=ll+1
                else if ( i==nnr*nnz ) then       !top right region
                    d(i-nnr)=aa
                    e(i-1)=aw
                    f(i)=-ap+aa
                    bi(i)=b(j+k)
                else if ( i==nj ) then         !left region
                    d(i-nnr)=aa
                    f(i)=-ap
                    g(i+1)=ae
                    h(i+nnr)=aa
                    bi(i)=b(j+k)
                    nj=nj+nnr
                    k=k+1
                    ll=ll+1
                else if ( i==nk ) then         !right region
                    d(i-nnr)=aa
                    e(i-1)=aw
                    f(i)=-ap
                    h(i+nnr)=aa
                    bi(i)=b(j+k)
                    nk=nk+nnr
                    j=j+nr
                    k=1
                    ll=2
                else if ( i < (nnr*nnz)-nnr ) then    !Interior region
                    d(i-nnr)=aa
                    e(i-1)=aw
                    f(i)=-ap
                    g(i+1)=ae
                    h(i+nnr)=aa
                    bi(i)=b(j+k)
                    k=k+1
                    ll=ll+1
                else                          !top region
                    d(i-nnr)=aa
                    e(i-1)=aw
                    f(i)=-ap+aa
                    g(i+1)=ae
                    bi(i)=b(j+k)
                    k=k+1
                    ll=ll+1
                end if  
            end do
            !For computational efficiency, the, A penta-diagonal matrix is split to 5 vectors
            call Bi_CGSTAB_P(d,e,f,g,h,bi,Vri) !Solves the linear system at each time-step
            j=nr+1
            k=nr+nr
            ll=1
            do i = 1,nr*nz !This loop adds the Bi_CGSTAB solved interior to the boundary radial velocity
                if ( i>nr+1 .and. i/=j .and. i/=k) then 
                    Vr(i)=Vri(ll)
                    ll=ll+1
                else if ( i==nr*nz-(nr-1) ) then !interior region
                    ll=nnr*nnz-(nnr-1)
                else if ( i==j ) then         !left region
                    j=j+nr
                else if ( i==k ) then         !right region
                    k=k+nr
                end if
            end do
            !print "(a5,25f10.7)",'Vr=', Vr(1:nr*nz)

            !The below loop solves the axial velocity at each grid point
            j=nr+1
            k=nr+nr
            ll=1
            do i = 1,nr*nz 
                aa=dz/(2*dr)
                ap=dz/ri(ll)
                k1=CLVZ-CLVZ*ri(ll)
                !The 9 if-statements below are for the 9 regions of our CV formulation
                if ( i==1 ) then !bottom left region
                    Vz(i)=k1 !inflow condition
                    ll=ll+1
                else if ( i < nr ) then        !bottom region
                    Vz(i)=k1
                    ll=ll+1
                else if ( i==nr ) then         !bottom right region
                    Vz(i)=0
                    ll=1
                else if ( i==nr*nz-(nr-1) ) then !top left region
                    Vz(i)=Vz(i-nr) !symmetry condition
                    ll=ll+1
                else if ( i==nr*nz ) then       !top right region
                    Vz(i)=0 !wall conditions
                else if ( i==j ) then         !left region
                    Vz(i)=Vz(i-nr) !symmetry conditions
                    j=j+nr
                    ll=ll+1
                else if ( i==k ) then         !right region
                    Vz(i)=0 !wall conditions
                    k=k+nr
                    ll=1
                else if ( i < nr*nz-nr ) then    !Interior region
                    Vz(i)=Vz(i-nr)-ap*Vr(i)-aa*Vr(i+1)+aa*Vr(i-1)
                    ll=ll+1
                else                          !top region
                    Vz(i)=Vz(i-nr) !outflow condition
                    ll=ll+1
                end if  
            end do
            !print "(a5,25f8.5)",'Vz=', Vz(1:nr*nz)
            k=0
            do j = 1, nz !this loop solves the mass flow-rate along each z-gridpoint
                mdot(j) = 0.d0
                do i = 1, nr
                    k=k+1
                    mdot(j) = mdot(j) + 2.d0*pi*rh*Vz(k)*ri(i)
                end do
            end do

            do i = 1,Nz
                som(i)=abs(mdot(i)-mdot(1))
            end do
            !print "(a5,25f8.3)",'som=', som(1:nz)
            !print "(a5,25f8.3)",'mdot=', mdot(1:nz)
            Mmdot=maxval(som)/mdot(1)
            Amdot=(sum(som)/nz)/mdot(1)

            !for printing results
            print "(i9,E12.3,E12.3)",ij, Mmdot, Amdot
            if (ij == 1) then
                Mmdot=1
                iii=101
            elseif (ij <= 100) then
                write(1,*) ij, Mmdot, Amdot
            elseif (ij == iii) then
                print *, achar(27)//"[2J" !clears the console
                iii=iii+5
            elseif (ij >= 3D6) then !set the maximum number of iteration
                Mmdot=tols
            end if

            !dt=dt+t !Adaptive time-step
            ij=ij+1
        end do

        close(1)
    end subroutine unsteady
!*************************************************************************************
    !*************************************************************************************
    subroutine Bi_CGSTAB_P(d,e,f,g,h,b,x)
        implicit none
        DOUBLE PRECISION, DIMENSION(nnr*nnz), INTENT(OUT) :: x
        DOUBLE PRECISION, dimension(nnr*nnz), INTENT(IN):: d, e, f, g, h, b
        DOUBLE PRECISION, DIMENSION(nnr*nnz) ::d1,e1,g1,h1,d2,e2,g2,h2
        DOUBLE PRECISION, DIMENSION(nnr*nnz) ::   p, r, r0, y, z, v, s, t, k, ks, kt
        DOUBLE PRECISION :: rho0, rho, w, alpha, beta, r_check, tol, nan
        INTEGER :: i
    !!This subroutine solves a penta-diagonal linear system using 
    !!The Preconditioned Bi_CGSTAB Algorithm by Van Der Vorst
    
        !Variable initialization
        d2=0; e2=0; g2=0; h2=0
        !since we are dealing with vectors, I have done the multiplication of two vectors, 
        !using this pattern
        !this multiples each vectors with x, for A*x
        d1=d(1:nnr*nnz)*x(1:nnr*nnz);e1=e(1:nnr*nnz)*x(1:nnr*nnz);g1=g(1:nnr*nnz)*x(1:nnr*nnz);h1=h(1:nnr*nnz)*x(1:nnr*nnz);
        !This circularlly shifts the vectors to the appropiate position, before addition
        d2(nnr+1:nnr*nnz)=d1(1:nnr*nnz-nnr);e2(2:nnr*nnz)=e1(1:nnr*nnz-1);g2(1:nnr*nnz-1)=g1(2:nnr*nnz)
        h2(1:nnr*nnz-nnr)=h1(nnr+1:nnr*nnz)
        !Hence A*x = d2+e2+f(1:nnr*nnz)*x(1:nnr*nnz)+g2+h2
        r0=b-(d2+e2+f(1:nnr*nnz)*x(1:nnr*nnz)+g2+h2)
        r=r0
        w=1; alpha=1; rho0=1; r_check=1
        v=0; p=0;  k=0
        tol=10E-8
        !For the inverse of K
        do i = 1, nnr*nnz
            k(i)=1/f(i) 
        end do
        i=0
        !main algorithm loop
        do while (r_check>tol)
            i=i+1
            rho=dot_product(r0,r)
            beta=(rho/rho0)*(alpha/w)
            rho0=dot_product(r0,r)
            p=r+beta*(p-w*v)
            y=k(1:nnr*nnz)*p(1:nnr*nnz)
            d1=d(1:nnr*nnz)*y(1:nnr*nnz);e1=e(1:nnr*nnz)*y(1:nnr*nnz);g1=g(1:nnr*nnz)*y(1:nnr*nnz);h1=h(1:nnr*nnz)*y(1:nnr*nnz);
            d2(nnr+1:nnr*nnz)=d1(1:nnr*nnz-nnr);e2(2:nnr*nnz)=e1(1:nnr*nnz-1);g2(1:nnr*nnz-1)=g1(2:nnr*nnz)
            h2(1:nnr*nnz-nnr)=h1(nnr+1:nnr*nnz)
            v=(d2+e2+f(1:nnr*nnz)*y(1:nnr*nnz)+g2+h2)
            alpha=rho/dot_product(r0,v)
            if (rho==0) then
                alpha=0
            end if
            s=r-alpha*v
            z=k(1:nnr*nnz)*s(1:nnr*nnz)
            d1=d(1:nnr*nnz)*z(1:nnr*nnz);e1=e(1:nnr*nnz)*z(1:nnr*nnz);g1=g(1:nnr*nnz)*z(1:nnr*nnz);h1=h(1:nnr*nnz)*z(1:nnr*nnz);
            d2(nnr+1:nnr*nnz)=d1(1:nnr*nnz-nnr);e2(2:nnr*nnz)=e1(1:nnr*nnz-1);g2(1:nnr*nnz-1)=g1(2:nnr*nnz)
            h2(1:nnr*nnz-nnr)=h1(nnr+1:nnr*nnz)
            t=(d2+e2+f(1:nnr*nnz)*z(1:nnr*nnz)+g2+h2)
            d1=d(1:nnr*nnz)*s(1:nnr*nnz);e1=e(1:nnr*nnz)*s(1:nnr*nnz);g1=g(1:nnr*nnz)*s(1:nnr*nnz);h1=h(1:nnr*nnz)*s(1:nnr*nnz);
            d2(nnr+1:nnr*nnz)=d1(1:nnr*nnz-nnr);e2(2:nnr*nnz)=e1(1:nnr*nnz-1);g2(1:nnr*nnz-1)=g1(2:nnr*nnz)
            h2(1:nnr*nnz-nnr)=h1(nnr+1:nnr*nnz)
            ks=(d2+e2+f(1:nnr*nnz)*s(1:nnr*nnz)+g2+h2)
            d1=d(1:nnr*nnz)*t(1:nnr*nnz);e1=e(1:nnr*nnz)*t(1:nnr*nnz);g1=g(1:nnr*nnz)*t(1:nnr*nnz);h1=h(1:nnr*nnz)*t(1:nnr*nnz);
            d2(nnr+1:nnr*nnz)=d1(1:nnr*nnz-nnr);e2(2:nnr*nnz)=e1(1:nnr*nnz-1);g2(1:nnr*nnz-1)=g1(2:nnr*nnz)
            h2(1:nnr*nnz-nnr)=h1(nnr+1:nnr*nnz)
            kt=(d2+e2+f(1:nnr*nnz)*t(1:nnr*nnz)+g2+h2)
            nan=dot_product(kt,ks)
            w=nan/dot_product(kt,kt)
            if (nan==0) then
                w=0
            end if
            x=x+alpha*y+w*z
            r=s-w*t 
            r_check=norm2(r)
        end do
    end subroutine Bi_CGSTAB_P
!*************************************************************************************
    
end program unsteady_Dfinite