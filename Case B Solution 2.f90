Module Parameters
  real(8),parameter::theta1=-90d0/180*3.1415926535897932384626433832795
  real(8),parameter::theta2=90d0/180*3.1415926535897932384626433832795
  real(8),parameter::nu=0.3
  real(8),parameter::pi=3.1415926535897932384626433832795
  integer,parameter::N=60
  integer,parameter::M11=2e4,M12=4e5
  integer,parameter::QQ=1e3
  real(8),parameter::p=1d0
  real(8),parameter::r=0.3
  real(8),parameter::epsilon=1d-20
  integer,parameter::T=45
  real(8),parameter::rho=0.5
  real(8),parameter::modulus=1
End Module Parameters

Subroutine calc
  use Parameters
  implicit none
  complex(8),dimension(2*N+1)::alpha,alpha1
  complex(8),dimension(2*N+1)::beta,beta1
  complex(8)::t1,t2,gamma,K0
  real(8)::kappa,lambda
  integer::k,l,i,j
  real(8),external::Factorial
  complex(8),external::h
  complex(8),dimension(N-1,N-1)::CoeffMat_n
  complex(8),dimension(N,N)::CoeffMat_p
  real(8),dimension(2*N+1)::ck11cos,ck11sin
  real(8),dimension(2*N+1)::ck12cos,ck12sin
  real(8),allocatable,dimension(:,:)::eta1,eta2
  complex(8),dimension(2*N+1)::ck11
  complex(8),dimension(2*N+1)::ck12
  real(8),dimension(N+1)::err_ck12_alpha
  real(8),dimension(N)::err_ck12_beta
  real(8),dimension(N+1)::err_ck11_alpha
  real(8),dimension(N)::err_ck11_beta
  real(8),allocatable,dimension(:)::tht
  complex(8),dimension(N-1)::dnq
  complex(8),dimension(N)::dpq
  complex(8),dimension(2*N+1)::dq
  complex(8),dimension(2*N+1)::d
  real(8)::p0,q0,pn1,qn1
  complex(8),dimension(N-1)::Constn
  complex(8),dimension(N)::Constp
  complex(8),dimension(2,2)::CoeffMat_01
  complex(8),dimension(2)::Const01
  external::ZGESVD,ZGESV
  complex(8),allocatable,dimension(:,:)::A
  complex(8),allocatable,dimension(:)::B
  real(8),allocatable,dimension(:)::S
  complex(8),allocatable,dimension(:,:)::U,VT
  integer,allocatable,dimension(:)::IPIV
  complex(8),allocatable,dimension(:)::WORK,RWORK
  integer::INFO
  real(8)::err
  integer::q
  complex(8),dimension(N+1)::Apk
  complex(8),dimension(N+1)::Bnk
  complex(8),dimension(2*N+1)::Ak
  complex(8),dimension(2*N+2)::Bk
  character(*),parameter::datapath=&
       '/home/luobinlin/researchpaper/paper6/data/'
  character(1)::var
  complex(8),dimension(T+1)::theta,sigma
  real(8),dimension(2*N+1)::Fk
  complex(8),dimension(T+1)::st
  real(8),dimension(T+1)::s2,s11,s22,t12
  complex(8),dimension(T+1)::str
  real(8),dimension(T+1)::s2r,s11r,s22r,t12r
  complex(8),dimension(T+1)::strho
  real(8),dimension(T+1)::s2rho,s11rho,s22rho,t12rho
  complex(8),dimension(T+1)::g
  real(8),dimension(T+1)::ux,uy
  complex(8),dimension(T+1)::gr
  real(8),dimension(T+1)::uxr,uyr
  complex(8),dimension(T+1)::grho
  real(8),dimension(T+1)::uxrho,uyrho

  t1= exp((0,1)*(theta1))
  t2= exp((0,1)*(theta2))
  kappa=3d0-4d0*nu
  lambda=log(kappa)/(2d0*pi)
  gamma=0.5+(0,1)*lambda
  K0=-(0,1)/2*exp(lambda/2*(theta1-theta2)-(0,1)/4*(theta1+theta2))

  alpha=0
  alpha1=0
  alpha(1)=-t1**(-gamma)*t2**(gamma-1)*1
  alpha(2)=-t1**(-gamma)*t2**(gamma-1)*(gamma*t1**(-1)-(gamma-1)*t2**(-1))
  do k=3,2*N+1
     do l=1,k-2
        alpha1(k)=alpha1(k)+h(-gamma,l)/Factorial(l)*h(gamma-1,k-1-l)/Factorial(k-1-l)&
             *t1**(-l)*t2**(-k+1+l)
     enddo
     alpha(k)=(-1)**(k-1)*h(-gamma,k-1)/Factorial(k-1)*t1**(-k+1)&
          +(-1)**(k-1)*h(gamma-1,k-1)/Factorial(k-1)*t2**(-k+1)&
          +(-1)**(k-1)*alpha1(k)
     alpha(k)=-t1**(-gamma)*t2**(gamma-1)*alpha(k)
  enddo

  beta=0
  beta1=0
  beta(1)=(1,0)
  beta(2)=gamma*t1-(gamma-1)*t2
  do k=3,2*N+1
     do l=1,k-2
        beta1(k)=beta1(k)&
             +h(-gamma,l)/Factorial(l)*h(gamma-1,k-l-1)/Factorial(k-l-1)*t1**(l)*t2**(k-l-1)
     enddo
     beta(k)=(-1)**(k-1)*h(-gamma,k-1)/Factorial(k-1)*t1**(k-1)&
          +(-1)**(k-1)*h(gamma-1,k-1)/Factorial(k-1)*t2**(k-1)&
          +(-1)**(k-1)*beta1(k)
  enddo

  ck12=0
  allocate(tht(M12),eta2(2*N+1,M12))
  tht=0
  do i=1,M12
     tht(i)=1d0*(i)/(M12+1)*(theta2-theta1)+theta1
  enddo
  ck12cos=0
  ck12sin=0
  eta2=0
  do k=-N,N
     do i=1,M12
        eta2(N+1+k,i)=1d0*(1d0*k+1d0/2)*tht(i)&
             +lambda*(log(sin((theta2-tht(i))/2))-log(sin((tht(i)-theta1)/2)))
     enddo
  enddo
  do k=-N,N
     do i=1,M12
        ck12cos(N+1+k)=ck12cos(N+1+k)+(sin((tht(i)-theta1)/2)*sin((theta2-tht(i))/2))**(-1d0/2)&
             *cos(eta2(N+1+k,i))*1d0/(M12+1)*(theta2-theta1)
     enddo
  enddo
  do k=-N,N
     do i=1,M12
        ck12sin(k+N+1)=ck12sin(N+1+k)+(sin((tht(i)-theta1)/2)*sin((theta2-tht(i))/2))**(-1d0/2)&
             *sin(eta2(N+1+k,i))*1d0/(M12+1)*(theta2-theta1)
     enddo
  enddo
  do k=-N,N
     ck12(N+1+k)=-exp(pi*lambda)*K0*ck12cos(N+1+k)-(0,1)*exp(pi*lambda)*K0*ck12sin(N+1+k)
  enddo
  open(1010,file=datapath//'/ck12/ck12_case_b.csv')
  do l=1,N
     write(1010,20002) l, ',' &
          , real(2d0*pi*(0,1)*kappa/(1d0+kappa)*alpha(l))/real(ck12(N+1-1-l+1)), ','&
          , aimag(2d0*pi*(0,1)*kappa/(1d0+kappa)*alpha(l))/aimag(ck12(N+1-1-l+1)), ','&
          , real(2d0*pi*(0,1)*kappa/(1d0+kappa)*beta(l))/real(ck12(N+1-1+l)), ','&
          , aimag(2d0*pi*(0,1)*kappa/(1d0+kappa)*beta(l))/aimag(ck12(N+1-1+l))
  enddo
  close(1010)
20002 format(I4,A,F10.6,A,F10.6,A,F10.6,A,F10.6)
  deallocate(tht,eta2)

  err_ck12_alpha=0
  err_ck12_beta=0
  do l=0,N-1
     err_ck12_alpha(l+1)=abs(alpha(l+1)+1d0*(1d0+kappa)/(2*pi*kappa*(0,1))*ck12(N+1-1-l))
  enddo
  do l=1,N
     err_ck12_beta(l)=abs(beta(l)-1d0*(1d0+kappa)/(2d0*pi*kappa*(0,1))*ck12(N+1-1+l))
  enddo
  write(*,*) 'Maximum error between ck12 and alpha =', maxval(err_ck12_alpha)
  write(*,*) 'Maximum error between ck12 and beta =', maxval(err_ck12_beta)

  allocate(tht(M11),eta1(2*N+1,M11))
  tht=0
  eta1=0
  do i=1,M11
     tht(i)=1d0*(i)/(M11+1)*((theta1-theta2+2d0*pi))+theta2
  enddo
  do k=-N,N
     do i=1,M11
        eta1(N+1+k,i)=1d0*(1d0*k+1d0/2)*tht(i)&
             +lambda*(log(sin((tht(i)-theta2)/2))-log(sin((tht(i)-theta1)/2)))
     enddo
  enddo
  ck11cos=0
  ck11sin=0
  do k=-N,N
     do i=1,M11
        ck11cos(k+N+1)=ck11cos(N+1+k)+(sin((tht(i)-theta1)/2)*sin((tht(i)-theta2)/2))**(-1d0/2)&
             *cos(eta1(N+1+k,i))*1d0/(M11+1)*(theta1-theta2+2d0*pi)
     enddo
  enddo
  do k=-N,N
     do i=1,M11
        ck11sin(k+N+1)=ck11sin(N+1+k)+(sin((tht(i)-theta1)/2)*sin((tht(i)-theta2)/2))**(-1d0/2)&
             *sin(eta1(N+1+k,i))*1d0/(M11+1)*(theta1-theta2+2d0*pi)
     enddo
  enddo
  ck11=0
  do k=-N,N
     ck11(N+1+k)=(0,1)*K0*ck11cos(N+1+k)-K0*ck11sin(N+1+k)
  enddo
  err_ck11_alpha=0
  err_ck11_beta=0
  do l=0,N-1
     err_ck11_alpha(l+1)=abs(kappa*alpha(l+1))/abs(ck11(N+1-1-l))
  enddo
  do l=1,N
     err_ck11_beta(l)=abs(beta(l))/abs(ck11(N+1-1+l))
  enddo
  write(*,*) 'Maximum error between ck11 and alpha =', maxval(err_ck11_alpha)-minval(err_ck11_alpha)
  write(*,*) 'Maximum error between ck11 and beta =', maxval(err_ck11_beta)-minval(err_ck11_beta)
  deallocate(tht,eta1)

  CoeffMat_n=0
  do i=1,N-1
     do j=1,N-1
        if (j >= i) then
           CoeffMat_n(i,j)=alpha(j-i+1)
        endif
     enddo
  enddo
  CoeffMat_p=0
  do i=1,N
     do j=1,N
        if (j >= i) then
           CoeffMat_p(i,j)=beta(j-i+1)
        endif
     enddo
  enddo
  CoeffMat_01=0
  CoeffMat_01(1,1)=alpha(1+0)
  CoeffMat_01(1,2)=-beta(1)
  CoeffMat_01(2,1)=ck11(N+1-1)
  CoeffMat_01(2,2)=ck11(N+1)

  allocate(A(N-1,N-1),S(N-1),U(N-1,N-1),VT(N-1,N-1),WORK(10*N),RWORK(5*N))
  A=CoeffMat_n
  S=0
  U=0
  VT=0
  WORK=0
  RWORK=0
  call ZGESVD('A','A',N-1,N-1,A,N-1,S,U,N-1,VT,N-1,WORK,10*N,RWORK,INFO)
  write(*,*) 'Condition number for Coefficient Matrix negative =', S(1)/S(N-1)
  deallocate(A,S,U,VT,WORK,RWORK)

  allocate(A(N,N),S(N),U(N,N),VT(N,N),WORK(10*(N)),RWORK(5*(N)))
  A=CoeffMat_p
  S=0
  U=0
  VT=0
  WORK=0
  RWORK=0
  call ZGESVD('A','A',N,N,A,N,S,U,N,VT,N,WORK,10*(N),RWORK,INFO)
  write(*,*) 'Condition number for Coefficient Matrix positive =', S(1)/S(N)
  deallocate(A,S,U,VT,WORK,RWORK)

  p0=0
  q0=0
  pn1=p
  qn1=0

  d=0
  Constn=0
  Constp=0
  Const01=0
  Constp(1)=r**2*(p0+(0,1)*q0)
  Constp(2)=1d0*2*(1-r**(-2))*r**3*(pn1-(0,1)*qn1)

  err=1d5
  q=0
  do while (err >= epsilon .and. q <= QQ)
     dnq=0
     dpq=0

     allocate(A(N-1,N-1),B(N-1),IPIV(N-1))
     A=CoeffMat_n
     B=Constn
     call ZGESV(N-1,1,A,N-1,IPIV,B,N-1,INFO)
     dnq=B
     deallocate(A,B,IPIV)

     allocate(A(N,N),B(N),IPIV(N))
     A=CoeffMat_p
     B=Constp
     call ZGESV(N,1,A,N,IPIV,B,N,INFO)
     dpq=B
     deallocate(A,B,IPIV)

     dq=0
     do k=1,N-1
        dq(k)=dnq(N-k)
     enddo
     do k=N+2,2*N+1
        dq(k)=dpq(k-N-1)
     enddo

     Const01=0
     do l=1,N-1
        Const01(1)=Const01(1)-alpha(l+1)*dq(N+1-1-l)
     enddo
     do l=2,N+1
        Const01(1)=Const01(1)+beta(l)*dq(N+1-1+l)
     enddo
     if (q==0) Const01(1)=Const01(1)+r*(pn1+(0,1)*qn1)
     do l=2,N
        Const01(2)=Const01(2)-ck11(N+1-l)*dq(N+1-l)
     enddo
     do l=1,N
        Const01(2)=Const01(2)-ck11(N+1+l)*dq(N+1+l)
     enddo
     allocate(A(2,2),B(2),IPIV(2))
     A=CoeffMat_01
     B=Const01
     call ZGESV(2,1,A,2,IPIV,B,2,INFO)
     dq(N)=B(1)
     dq(N+1)=B(2)
     deallocate(A,B,IPIV)

     Apk=0
     Bnk=0
     do k=0,N
        do l=0,N+k
           Apk(k+1)=Apk(k+1)+alpha(l+1)*dq(N+1+k-l)
        enddo
     enddo
     do k=1,N
        do l=1,N+k
           Bnk(k)=Bnk(k)+beta(l)*dq(N+1-k+l)
        enddo
     enddo

     Constn=0
     Constp=0
     do k=2,N
        Constn(k-1)=1d0*(k-1)*(1-r**(-2))*r**(2*k)*conjg(Apk(k+1))+r**(2*k-2)*Bnk(k)
     enddo
     Constp(1)=r**2*Apk(0+1)+(1-r**(-2))*r**2*conjg(Apk(0+1))
     Constp(2)=r**4*Apk(1+1)+2d0*(1-r**(-2))*r**2*conjg(Bnk(1))
     do k=2,N-1
        Constp(k+1)=r**(2*k+2)*Apk(k+1)+1d0*(k**2-1)*(1-r**(-2))**2*r**(2*k+2)*Apk(k+1)&
             +1d0*(k+1)*(1-r**(-2))*r**(2*k)*conjg(Bnk(k))
     enddo

     err=maxval(abs(dq))
     q=q+1
     d=d+dq
     ! write(*,*) 'Error of', q, 'th iteration =',err
  enddo
  ! d=-real(d)+(0,1)*aimag(d)
  write(*,*) 'Error of', q, 'th iteration =',err
  open(2001,file=datapath//'/dk/case_b_solution2.csv')
  do k=-N,N
     write(2001,20001) k, ',', real(d(N+1+k)), ',', aimag(d(N+1+k))
  enddo
  close(2001)
20001 format(I4,A,F10.6,A,F10.6)

  Ak=0
  Bk=0
  do k=-N,N
     do l=0,N+k
        Ak(k+N+1)=Ak(k+N+1)+alpha(l+1)*d(k-l+N+1)
     enddo
  enddo
  do k=-N-1,N-1
     do l=1,N-k
        Bk(k+N+2)=Bk(k+N+2)+beta(l)*d(k+l+N+1)
     enddo
  enddo
  Bk(N+N+2)=0

  ! write(*,*) 'Error of A(-1) =', abs(Ak(N)-p*r/(1d0+kappa))/abs(p*r/(1d0+kappa))
  ! write(*,*) 'Error of B(-1) =', abs(Bk(N+1)+kappa*p*r/(1d0+kappa))/abs(kappa*p*r/(1d0+kappa))

  do k=-N,N
     Fk(k+N+1)=sin((abs(k))*pi/(N))/((abs(k))*pi/(N))
  enddo
  Fk(N+1)=1d0

  do i=1,T+1
     theta(i)=1d0*(i-1)/T*2d0*pi
     sigma(i)=exp((0,1)*theta(i))
  enddo
  st=0
  do i=1,T+1
     do k=-N,N
        st(i)=st(i)+(Ak(k+N+1)-Bk(k+N+2))*sigma(i)**(k)*Fk(k+N+1)
     enddo
  enddo
  s11=real(st)
  t12=aimag(st)

  str=0
  do i=1,T+1
     do k=-N,N
        str(i)=str(i)+&
             (Ak(k+N+1)*r**(k)&
             +1d0*(k+1)*conjg(Ak(-k+N+1))*(1-r**(-2))*r**(-k)&
             -Bk(k+N+2)*r**(-k-2))*sigma(i)**(k)*Fk(k+N+1)
     enddo
  enddo
  s11r=real(str)
  t12r=aimag(str)

  strho=0
  do i=1,T+1
     do k=-N,N
        strho(i)=strho(i)+&
             (Ak(k+N+1)*rho**(k)&
             +1d0*(k+1)*conjg(Ak(-k+N+1))*(1-rho**(-2))*rho**(-k)&
             -Bk(k+N+2)*rho**(-k-2))*sigma(i)**(k)*Fk(k+N+1)
     enddo
  enddo
  s11rho=real(strho)
  t12rho=aimag(strho)

  s2=0
  do i=1,T+1
     do k=-N,N
        s2(i)=s2(i)+4d0*real(Ak(k+N+1)*sigma(i)**k*Fk(k+N+1))
     enddo
  enddo
  s22=s2-s11

  s2r=0
  do i=1,T+1
     do k=-N,N
        s2r(i)=s2r(i)+4d0*real(Ak(k+N+1)*r**k*sigma(i)**k*Fk(k+N+1))
     enddo
  enddo
  s22r=s2r-s11r

  s2rho=0
  do i=1,T+1
     do k=-N,N
        s2rho(i)=s2rho(i)+4d0*real(Ak(k+N+1)*rho**k*sigma(i)**k*Fk(k+N+1))
     enddo
  enddo
  s22rho=s2rho-s11rho

  open(1001,file=datapath//'/case_b/stress_case_b2.csv')
  do i=1,T+1
     write(1001,10001) 1d0*(i-1)/T*360, ',', &
          -s11(i), ',', -s11r(i), ',', -s11rho(i), ',', &
          -s22(i), ',', -s22r(i), ',', -s22rho(i), ',', &
          -t12(i), ',', -t12r(i), ',', -t12rho(i)
  enddo
  close(1001)

  g=0
  do i=1,T+1
     do k=0,N
        g(i)=g(i)+1d0*(kappa*Ak(k+N+1)/(k+1)&
             -Bk(k+N+2)/(-k-1))*sigma(i)**(k+1)*Fk(N+1+k)
     enddo
     g(i)=g(i)-conjg(Ak(N+1+1))*Fk(N+1+1)
     do k=2,N
        g(i)=g(i)+1d0*(kappa*Ak(N+1-k)/(-k+1)&
             -Bk(N+2-k)/(k-1))*sigma(i)**(-k+1)*Fk(N+1+k)
     enddo
  enddo
  ux=real(g)/(2d0*modulus)
  uy=aimag(g)/(2d0*modulus)

  gr=0
  do i=1,T+1
     do k=0,N
        gr(i)=gr(i)+1d0*(kappa*Ak(N+1+k)*r**(k+1)/(k+1)&
             -conjg(Ak(N+1-k))*r**(-k+1)*(1-r**(-2)) &
             -Bk(k+N+2)*r**(-k-1)/(-k-1))*sigma(i)**(k+1)*Fk(N+1+k)
     enddo
     gr(i)=gr(i)+(kappa*Ak(N)*log(r)-conjg(Ak(N+1+1))*r**2-Bk(N+1)*log(r))*Fk(N+1+1)
     do k=2,N
        gr(i)=gr(i)+(kappa*Ak(N+1-k)*r**(-k+1)/(-k+1)&
             -conjg(Ak(N+1+k))*r**(k+1)*(1-r**(-2)) &
             -Bk(N+2-k)*r**(k-1)/(k-1))*sigma(i)**(-k+1)*Fk(N+1+k)
     enddo
  enddo
  uxr=real(gr)/(2d0*modulus)
  uyr=aimag(gr)/(2d0*modulus)

  grho=0
  do i=1,T+1
     do k=0,N
        grho(i)=grho(i)+(kappa*Ak(k+N+1)*rho**(k+1)/(k+1)&
             -conjg(Ak(N+1-k))*rho**(-k+1)*(1-rho**(-2)) &
             -Bk(N+2+k)*rho**(-k-1)/(-k-1))*sigma(i)**(k+1)*Fk(N+1+k)
     enddo
     grho(i)=grho(i)+(kappa*Ak(N)*log(rho)-conjg(Ak(N+1+1))*rho**2-Bk(N+1)*log(rho))*Fk(N+1+1)
     do k=2,N
        grho(i)=grho(i)+(kappa*Ak(N+1-k)*rho**(-k+1)/(-k+1)&
             -conjg(Ak(N+1+k))*rho**(k+1)*(1-rho**(-2)) &
             -Bk(N+2-k)*rho**(k-1)/(k-1))*sigma(i)**(-k+1)*Fk(N+1+k)
     enddo
  enddo
  uxrho=real(grho)/(2d0*modulus)
  uyrho=aimag(grho)/(2d0*modulus)

  open(1002,file=datapath//'/case_b/displacement_case_b2.csv')
  do i=1,T+1
     write(1002,10002) 1d0*(i-1)/T*360, ',', &
          -(ux(i)-ux(T))*1e1, ',', -(uxr(i)-ux(T))*1e1, ',', -(uxrho(i)-ux(T))*1e1, ',', &
          -(uy(i)-uy(T))*1e1, ',', -(uyr(i)-uy(T))*1e1, ',', -(uyrho(i)-uy(T))*1e1
  enddo
  close(1002)
  ! do i=1,T
  !    write(*,*) g(i),gr(i),grho(i)
  ! enddo
  write(*,*) 'Horizontal rigid-body displacement of Case B Sol 2 =', ux(T)*1e1
  write(*,*) 'Vertical rigid-body displacement of Case B Sol 2 =', uy(T)*1e1

10001 format(F12.7,A,&
           F12.7,A,F12.7,A,F12.7,A,&
           F12.7,A,F12.7,A,F12.7,A,&
           F12.7,A,F12.7,A,F12.7)
10002 format(F12.7,A,&
           F12.7,A,F22.7,A,F12.7,A,&
           F12.7,A,F22.7,A,F12.7)
End Subroutine calc

Function h(in,k)
  implicit none
  complex(8)::h
  complex(8),intent(in)::in
  integer,intent(in)::k
  integer::i

  h=in
  do i=1,k-1
     h=h*(in-i)
  enddo
End Function h

Function Factorial(k)
  implicit none
  real(8)::Factorial
  integer,intent(in)::k
  integer::i

  Factorial=1d0
  do i=1,k
     Factorial=Factorial*i
  enddo
End Function Factorial

Program Main
  implicit none
  call calc
End Program Main
