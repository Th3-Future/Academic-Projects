!****************Begin Header***************************************************
!This program was written by Godswill Ezeorah, Student Number: 501012886 on June 10, 2021.
!This program solves a linear/non-linear equation using finite volume method
!and was written as a solution to AE8112 PS2 q1a
!****************************End Header*****************************************
program newt_multiv
    implicit none
    DOUBLE PRECISION, dimension(:), ALLOCATABLE :: T_bar, f_T, Ti, ri
    DOUBLE PRECISION, dimension(:,:), ALLOCATABLE :: Jac
    DOUBLE PRECISION, PARAMETER :: Pi = 4*atan(1.0), del=8.4*10.0**(-8), r=0.10
    DOUBLE PRECISION :: tol, norm, Pin, Pout, Pint, HL_error, Teg, start, finish
    INTEGER :: n, i1
    call cpu_time(start) !gets the start time, for timing purpose
    !Initialization of variables
    tol=10.0**(-9) !The given tolerance
    Teg=800        !Temperature at the edge of the disc
    open(1, file = 'PS2_Q1a.csv', status = 'unknown')
    !Loop through for n=100, 200 and 300 control volumes
    do n = 100, 300, 100
        ALLOCATE(T_bar(n),f_T(n),jac(n,n),Ti(n),ri(n))  
        norm=1 
        i1=1
        !Initial geuss
        T_bar=Teg
        do while (norm>tol)
            call solve_fT(T_bar,F_T, ri)
            call solve_jx(jac, T_bar)
            call solve_gauss(jac,f_T,Ti)
            T_bar=T_bar+Ti
            norm=norm2(Ti)/n
            i1=i1+1
        end do
        Pint=pi*0.001*(Teg-T_bar(n))*4*1000 !Another approach to calculating Power-in
        call solve_P(T_bar,Pin,Pout)
        HL_error=100*abs(Pin-Pout)/Pin  !calculates % heat loss error
        !Printing our results
        print *, "no. of Control Volume", n
        print *, "no. Iteration before Convergence", i1
        write(*,2) "T(r) = ", T_bar
        write(*,3) "Power in = ", Pin,"W"
        write(*,3) "Pin another_method = ", Pint,"W"
        write(*,3) "Power out = ", Pout,"W"
        write(*,3) "% Heat Loss Error = ", HL_error
        if (n<300) then
            DEALLOCATE(T_bar,f_T,jac,Ti,ri)
        end if 
    end do
    do i1 = 1, 300
        write(1,*) ri(i1), T_bar(i1) !Writing the results to a .txt file
    end do
    close(1)
    call cpu_time(finish) !gets the end time
    write(*,1) "Program Execution Time = ", (finish-start)/60,"min"
    1 format(a40,f10.3,a3)
    2 format(a7,300f10.3)
    3 format(a40,f10.3,a2)

    contains
!*************************************************************************************
    subroutine solve_fT(Tr, fT, rs)
        DOUBLE PRECISION, dimension(n), INTENT(OUT) ::  fT, rs
        DOUBLE PRECISION, dimension(n) :: Tr
        DOUBLE PRECISION :: aa, ab, TW, TP, TE, rr, dr
        INTEGER :: i
    !!This subroutine creates and solve the discretized function
        
        !Initialization of variables
        dr=r/n
        rr=r
        aa=1000/dr    !for equispaced grid
        ab=(11.34*10.0**(-5))/16
        do i = 1,n !Loops through our control volumes to solve the function vector
            TP=Tr(i)
            if ( i==1 ) then
                TE=Tr(i+1)
                fT(1)=-aa*(TE-TP) + 2*aa*(TP-Teg) - ab*(rr-dr)*(TE+TP)**4 + ab*16*rr*Teg**4
            else if ( i<n ) then
                TW=Tr(i-1)
                TE=Tr(i+1)
                fT(i)=-aa*(TE-TP) + aa*(TP-TW) - ab*(rr-dr)*(TE+TP)**4 + ab*rr*(TP+TW)**4
            else
                TW=Tr(i-1)
                fT(n)= aa*(TP-TW) - ab*16*(rr-dr)*(TP)**4 + ab*rr*(TP+TW)**4
            end if
            rs(i)=rr-dr/2
            rr=rr-dr
        end do
        fT=-fT !Using Newton's methos, this has to be negative
    end subroutine solve_fT
!*************************************************************************************
!*************************************************************************************
    subroutine solve_Jx(Jx, T_b)
        DOUBLE PRECISION, dimension(n,n) ::  Jx
        DOUBLE PRECISION, dimension(n) :: f_per, T_b, rs
        DOUBLE PRECISION :: per
        integer  :: j, k
    !!This subroutine creates and solve the Jacobian of our function
    !!Using numerically evaluated Jacobian

        do j=1,n
            per=del*T_b(j)+del
            T_b(j)=T_b(j)+per
            call solve_fT(T_b, f_per, rs)
            T_b(j)=T_b(j)-per
            do k=1,n
                Jx(k,j)=(-f_per(k)+f_T(k))/per
            end do
        end do
    end subroutine solve_Jx
!*************************************************************************************
!*************************************************************************************
    subroutine solve_P(Tr, P1, P2)
        DOUBLE PRECISION, INTENT(OUT) :: P1, P2
        DOUBLE PRECISION, dimension(n) ::  qin, qout
        DOUBLE PRECISION, dimension(n) :: Tr
        DOUBLE PRECISION :: aa, ac, ab, TW, TP, rr, dr
        INTEGER :: i
    !!This subroutine solves the Power-in and Power-out, using our formulation

        !Variable initialization
        dr=r/n
        rr=r
        aa=0.001*pi
        ac=pi*11.34*10.0**(-8)   
        ab=aa*4*1000
        TW=Teg   !Temperature maintained at the edge
        do i = 1,n
            TP=Tr(i)
            if ( i==1 ) then
                qin(1)=ab*(rr**2-(rr-dr)**2)*(TW-Tp)/(rr**2-(rr-dr)**2)
                qout(1)=ac*(rr**2-(rr-dr)**2)*TP**4
            else if ( i<n ) then
                TW=Tr(i-1)
                qin(i)=ab*(rr**2-(rr-dr)**2)*(TW-Tp)/(rr**2-(rr-dr)**2)
                qout(i)=ac*(rr**2-(rr-dr)**2)*TP**4
            else
                TW=Tr(i-1)
                qin(n)=ab*(rr**2-(rr-dr)**2)*(TW-Tp)/(rr**2-(rr-dr)**2)
                qout(n)=ac*(rr**2-(rr-dr)**2)*TP**4
            end if
            rr=rr-dr
        end do
        P1=sum(qin)
        P2=2*sum(qout) !Multiplied by 2, because radiation occurs from the top and bottom.
    end subroutine solve_P
!*************************************************************************************
!*************************************************************************************
    subroutine solve_gauss(Ag1, b, X)
        implicit none
        DOUBLE PRECISION, dimension(n), INTENT(OUT)   :: X
        DOUBLE PRECISION, dimension(n,n)              :: Ag1
        DOUBLE PRECISION, dimension(n,n+1)            :: A
        DOUBLE PRECISION, dimension(n+1)              :: prod, new, switch1, switch2
        DOUBLE PRECISION, dimension(n)                :: b
        INTEGER, dimension(1)                     :: ros
        DOUBLE PRECISION                              :: ratio
        INTEGER                                   :: steps, k1, m1, t1, i, j
        !! This subroutine solve a linear problem using Guassian Elimination method (with partial pivoting)
        
        !Initialization of variables
        X=0
        t1=1
        !rewriting the given matrices in the augumented matrix form
        A(1:n,1:n)=Ag1
        A(1:n,n+1)=b
        !!Foreward elimination
        steps=n-1                      !Total number of forward elimination steps
        do i=1,steps                   !for new pivot point (diagonally, starting from 1,1)
            ros=maxloc(abs(A(i:n,i)))
                
            if ( ros(1)==1 ) then       !if the 1st element of the 1st row is the largest (i.e. at pivot)
                do j=i,steps           !for consecutive rows below the pivot
                    ratio=A(j+1,i)/A(i,i)
                    prod=A(i,:)*ratio
                    new=A(j+1,:)-prod
                    A(j+1,:)=new
                end do
                else                 !else switch the 1st row with the with the row with the largest column element
                switch1=A(i,:)   
                switch2=A(ros(1),:)
                A(i,:)=switch2
                A(ros(1),:)=switch1
                do j=i,steps              !Steps
                    ratio=A(j+1,i)/A(i,i)
                    prod=A(i,:)*ratio
                    new=A(j+1,:)-prod
                    A(j+1,:)=new
                end do
                end if 
        end do
      
        !!Since the upper diagonal is obtained from the above, the X varibles can be solved as:
        ! x8=A(n,n+1)/A(n,n)
        ! x7=A(n,n+1)/A(n-1,n-1)-x8*A(n-1,n)/A(n-1,n-1)
        ! x6=A(n,n+1)/A(n-2,n-2)-x8*A(n-2,n)/A(n-2,n-2)-x7*A(n-2,n-1)/A(n-2,n-2)
        ! x5=A(n,n+1)/A(n-3,n-3)-x8*A(n-3,n)/A(n-3,n-3)-x7*A(n-3,n-1)/A(n-3,n-3)-x6*A(n-3,n-2)/A(n-3,n-3)
        !...
        !using the above pattern we can create the backward substition, using a loop and a recuresive procedures
        !as shown below
      
      !! this function account for the vertical change in the above pattern
        do k1 = steps, 0, -1
            m1=n-k1
            if ( k1==0 ) then
                 X(n)=A(n,n+1)/A(n,n)
             else
                 X(n)=A(n,n+1)/A(n,n)
                 X(k1)=A(k1,n+1)/A(k1,k1)-r_fun(A,X,k1,m1,t1)                
            end if
        end do  
        end subroutine solve_gauss
      !! this function accounts for the horinzontal change in the pattern
              recursive function r_fun(A1,Xr,k,m,t) result(fr)
                  DOUBLE PRECISION, dimension(n,n+1) :: A1
                  DOUBLE PRECISION, dimension(n) :: xr
                  DOUBLE PRECISION :: fr
                  INTEGER :: m, k, t
                  if ( m==0 ) then
                      fr=Xr(n-(t-1))*(A1(k,n-(t-1))/A1(k,k))
                  else
                      fr=Xr(n-(t-1))*(A1(k,n-(t-1))/A1(k,k))+r_fun(A1,Xr,k,m-1,t+1)
                  end if
              end function r_fun
!*************************************************************************************
end program newt_multiv