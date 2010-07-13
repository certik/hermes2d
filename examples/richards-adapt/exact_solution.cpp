
// funkce pro vypocet hodnoty na okrajove podmince, mimo hodnotu jeste spocte derivaci funkce na teto hrane vuci x_1, ale to je asi zbytecny
double boundary_value(double x, double& dhdx) {
//  a...intent(in)..materialovy parametr
//  H_R.intent(in)..saci tlak kdyz je material "very dry"
//  x.. intent(in)...x_1-slozka bodu na hranici \Gamma_3
//  dh..intent(out)..derivace okrajove podminky...
//  h...intent(out)..hodnota okrajove podminky na hranici \Gamma_3
  dhdx = (1-exp(ALPHA*H_R)*M_PI*cos(M_PI*x/ALPHA))/(ALPHA*ALPHA*(exp(A*H_R)+(1-exp(ALPHA*H_R))*sin(M_PI*x/ALPHA))); //tohle asi nepotrebujes, derivace op podle x_1 na hrane
  return 1/ALPHA*log(exp(ALPHA*H_R) + (1-exp(ALPHA*H_R))*sin(M_PI*x/ALPHA)); //hodnota bodu na hrane s \Gamma_3
}


// funkce pro vypocet funkcni hondoty reseni a derivace reseni na domene
  double exact_sol(double x, double z, double& dhdx, double& dhdz) {
// !--------i/o variables-------------------------
//       !> coordinates of the desired point
//       real(kind=rkind), dimension(:), intent(in) :: x,z
//       !> simulation time to plot
//       real(kind=rkind), intent(in)               :: TIME
//       !> material parameter
//       real(kind=rkind), dimension(:), intent(in) :: ALPHA
//       !> "a very dry" media pressure head
//       real(kind=rkind, intent(in)                :: H_R
//       !> saturated water content(porosity) and residual water content
//       real(kind=rkind), intent(in)               :: THETA_R, THETA_S
//       !> saturated hydraulic conductivity
//       real(kind=rkind), intent(in)               :: K_S
//       !> domain length(x_3) + width(x_1)   --- tohle jsem vlastne nezkousel kdyz to bude jiny tvar nez ctverec, muze se stat, ze to bude treba prohodit
//       real(kind=rkind), intent(in)               :: L, A 
//       !> pressure head at the desired point
//       real(kind=rkind), intent(out)              :: h
//       !> solution derivations
//       real(kind=rkind), intent(out)              :: dhdx, dhdz


// !--------local variables--------------------
   
  
      
// pro kontrolu seznam lokalnich promennych z originalniho kodu
// 
//       real(kind=rkind) :: lambda
//       real(kind=rkind) :: c
//       real(kind=rkind) :: gamma
//       real(kind=rkind) :: phi
//       real(kind=rkind) :: hbar
//       real(kind=rkind) :: beta
//       real(kind=rkind) :: ho
//       real(kind=rkind) :: H_R
//       real(kind=rkind) :: suma
//       real(kind=rkind) :: hss
//       real(kind=rkind) :: tmp
//       real(kind=rkind) :: ALPHA
//       real(kind=rkind) :: A
//       real(kind=rkind) :: L, x,z,pi, a_const, f_xz, g_xz, k_xz, dfdx, dfdz, dgdx, dgdz, dkdx, dkdz, suma2
//       integer :: i


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !-----end of declarations--------------------------------------------------
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


// !----formula parameters
      double reps = 1e-9 ; //presnost realnych cisel

      double ho = 1 - exp(ALPHA*H_R) ;

      double beta = sqrt(ALPHA*ALPHA/4 + (M_PI/A)*(M_PI/A)) ;

      double c = ALPHA*(THETA_S - THETA_R)/K_S ;

      double hss = ho*sin(M_PI*x/A)*exp(ALPHA/2*(L-z))*sinh(beta*z)/sinh(beta*L) ;

      double sum = 0 ;


// !----fourier series for solution and for solution x_1 derivate
      for (int i=1; i >= 0; i++) {
        double lambda = i*M_PI/L ;
        double gamma = 1/c*(beta*beta + lambda*lambda) ;
        double tmp = pow(-1,i)*lambda/gamma*sin(lambda*z)*exp(-gamma*TIME) ;
        sum += tmp ;
        if (fabs(1/c*lambda/gamma * exp(-gamma*TIME)) < reps*reps) break;
      }  


     double phi = 2*ho/(L*c)*sin(M_PI*x/A)*exp(ALPHA/2*(L-z))*sum ;


     double hbar = phi + hss ;

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!! vysledna hodnota funkce h!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
     double h = 1/ALPHA*log(exp(ALPHA*H_R)+hbar) ; // !!!!! function h !!!!
//      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
//      !!!!!!vypocet derivaci!!!!!!!!!

     double a_const = exp(ALPHA*H_R) ;

     double f_xz = 2*ho/(L*c)*sin(M_PI*x/A)*exp(ALPHA/2*(L-z)) ; 

     double g_xz = sum ;

     double k_xz = ho*sin(M_PI*x/A)*exp(ALPHA/2*(L-z))*sinh(beta*z)/sinh(beta*L) ;
     
     double dfdx = 2*exp(1/2*ALPHA*(L-z))*ho*M_PI*cos(M_PI*x/A)/(A*c*L) ;

     double dfdz = -ALPHA*exp(1/2*ALPHA*(L-z))*ho*sin(M_PI*x/A)/(c*L) ;

     double sum2 = 0 ;
     
     for (int i=1; i >= 0; i++) {
       double lambda = i*M_PI/L;
       double gamma = 1/c*(beta*beta + lambda*lambda);
       double tmp  = pow(-1,i)*i*exp(-gamma*TIME)*lambda*lambda*cos(lambda*z)/gamma;
       sum2 += tmp;
       if (fabs(lambda*lambda*exp(-gamma*TIME)/gamma) < reps*reps) break;
     }
       
       
    double dgdz = sum2 ;
    
    double dgdx = 0 ;
    
    double dkdx = exp(1/2*ALPHA*(L-z))*ho*M_PI*cos(M_PI*x/A)*sinh(beta*z)/(A*sinh(beta*L)) ;

    double dkdz = beta*exp(1/2*ALPHA*(L-z))*ho*cosh(beta*z)*sin(M_PI*x/A)/sinh(beta*L)- ALPHA*exp(1/2*ALPHA*(L-z))*ho*sin(M_PI*x/A)*sinh(beta*z)/(2*sinh(beta*L)) ;
    
//  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
// !!!!!!!!!vysledne hodnoty derivaci!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
    dhdx = 1/ALPHA*1/(a_const + f_xz*g_xz +k_xz)*(dfdx*g_xz + dgdx*f_xz + dkdx) ;
    
    dhdz = 1/ALPHA*1/(a_const + f_xz*g_xz +k_xz)*(dfdz*g_xz + dgdz*f_xz + dkdz) ;
    
    //  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
// !!!!!!!!!vysledne hodnoty derivaci!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111

    return h;

  }
