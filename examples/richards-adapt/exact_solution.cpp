
// funkce pro vypocet hodnoty na okrajove podmince, mimo hodnotu jeste spocte derivaci funkce na teto hrane vuci x_1, ale to je asi zbytecny
double boundary_value(double h_r, double x, double& dh, double& h) {
//  a...intent(in)..materialovy parametr
//  h_r.intent(in)..saci tlak kdyz je material "very dry"
//  x.. intent(in)...x_1-slozka bodu na hranici \Gamma_3
//  dh..intent(out)..derivace okrajove podminky...
//  h...intent(out)..hodnota okrajove podminky na hranici \Gamma_3
  dh = (1-exp(ALPHA*h_r)*M_PI*cos(M_PI*x/ALPHA))/(ALPHA*ALPHA*(exp(a*h_r)+(1-exp(ALPHA*h_r))*sin(M_PI*x/ALPHA))); //tohle asi nepotrebujes, derivace op podle x_1 na hrane
  h = 1/ALPHA*log(exp(ALPHA*h_r)+(1-exp(ALPHA*h_r))*sin(M_PI*x/ALPHA)); //hodnota bodu na hrane s \Gamma_3
}


// funkce pro vypocet funkcni hondoty reseni a derivace reseni na domene
  double exact_sol(double x, double z, double hr, double a, double L, double& dhdx, double& dhdz) {
// !--------i/o variables-------------------------
//       !> coordinates of the desired point
//       real(kind=rkind), dimension(:), intent(in) :: x,z
//       !> simulation time to plot
//       real(kind=rkind), intent(in)               :: TIME
//       !> material parameter
//       real(kind=rkind), dimension(:), intent(in) :: ALPHA
//       !> "a very dry" media pressure head
//       real(kind=rkind, intent(in)                :: hr
//       !> saturated water content(porosity) and residual water content
//       real(kind=rkind), intent(in)               :: THETA_R, THETA_S
//       !> saturated hydraulic conductivity
//       real(kind=rkind), intent(in)               :: K_S
//       !> domain length(x_3) + width(x_1)   --- tohle jsem vlastne nezkousel kdyz to bude jiny tvar nez ctverec, muze se stat, ze to bude treba prohodit
//       real(kind=rkind), intent(in)               :: L, a 
//       !> pressure head at the desired point
//       real(kind=rkind), intent(out)              :: h
//       !> solution derivations
//       real(kind=rkind), intent(out)              :: dhdx, dhdz


// !--------local variables--------------------
      double lambda ;
      double c ;
      double gamma ;
      double phi ;
      double hbar ;
      double beta ;
      double ho ;
      double suma ;
      double suma2 ;
      double hss ;
      double tmp ;
      double reps ;
      double a_const ;
      double f_xz ;
      double g_xz ;
      double k_xz ;
      double dfdx ;
      double dfdz ;
      double dgdx ;
      double dgdz ;
      double dkdx ;
      double dkdz ;
      int i;
  
      
// pro kontrolu seznam lokalnich promennych z originalniho kodu
// 
//       real(kind=rkind) :: lambda
//       real(kind=rkind) :: c
//       real(kind=rkind) :: gamma
//       real(kind=rkind) :: phi
//       real(kind=rkind) :: hbar
//       real(kind=rkind) :: beta
//       real(kind=rkind) :: ho
//       real(kind=rkind) :: hr
//       real(kind=rkind) :: suma
//       real(kind=rkind) :: hss
//       real(kind=rkind) :: tmp
//       real(kind=rkind) :: ALPHA
//       real(kind=rkind) :: a
//       real(kind=rkind) :: L, x,z,pi, a_const, f_xz, g_xz, k_xz, dfdx, dfdz, dgdx, dgdz, dkdx, dkdz, suma2
//       integer :: i


// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !-----end of declarations--------------------------------------------------
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


// !----formula parameters
      reps = 1e-9 ; //presnost realnych cisel

      ho = 1-exp(ALPHA*hr) ;

      beta = sqrt(ALPHA**2/4 + (M_PI/a)**2) ;

      c = ALPHA*(parameters(5)-parameters(4))/parameters(6) ;

      hss = ho*sin(M_PI*x/a)*exp(ALPHA/2*(L-z))*sinh(beta*z)/sinh(beta*L) ;

      suma = 0 ;

      i = 0 ;


// !----fourier series for solution and for solution x_1 derivate
      do {
        i = i+1 ;
        tmp = suma ;
        lambda = i*M_PI/L ;
        gamma = 1/c*(beta*beta + lambda*lambda) ;
        tmp = pow(-1,i)*lambda/gamma*sin(lambda*z)*exp(-gamma*TIME) ;
        suma = suma + tmp ;
      } while (abs(1/c*lambda/gamma * exp(-gamma*TIME)) > reps*reps);  


     phi = 2*ho/(L*c)*sin(M_PI*x/a)*exp(ALPHA/2*(L-z))*suma ;


     hbar = phi + hss ;

// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!! vysledna hodnota funkce h!!!!!!!!!!!!!!!!!!!!!!!!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
     h = 1/ALPHA*log(exp(ALPHA*hr)+hbar) ; // !!!!! function h !!!!
//      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     
     
//      !!!!!!vypocet derivaci!!!!!!!!!

     a_const = exp(ALPHA*hr) ;

     f_xz = 2*ho/(L*c)*sin(M_PI*x/a)*exp(ALPHA/2*(L-z)) ; 

     g_xz = suma ;

     k_xz = ho*sin(M_PI*x/a)*exp(ALPHA/2*(L-z))*sinh(beta*z)/sinh(beta*L) ;
     
     dfdx = 2*exp(1/2*ALPHA*(L-z))*ho*M_PI*cos(M_PI*x/a)/(a*c*L) ;

     dfdz = -ALPHA*exp(1/2*ALPHA*(L-z))*ho*sin(M_PI*x/a)/(c*L) ;

     i = 0 ;
     suma2 = 0 ;
     
     do{
       i = i+1;
       lambda = i*M_PI/L;
       gamma = 1/c*(beta*beta + lambda*lambda);
       tmp  = pow(-1,i)*i*exp(-gamma*TIME)*lambda*lambda*cos(lambda*z)/gamma;
       suma2 = suma2 + tmp;
     } while (abs(lambda*lambda*exp(-gamma*TIME)/gamma) > reps*reps);
       
       
    dgdz = suma2 ;
    
    dgdx = 0 ;
    
    dkdx = exp(1/2*ALPHA*(L-z))*ho*M_PI*cos(M_PI*x/a)*sinh(beta*z)/(a*sinh(beta*L)) ;

    dkdz = beta*exp(1/2*ALPHA*(L-z))*ho*cosh(beta*z)*sin(M_PI*x/a)/sinh(beta*L)- ALPHA*exp(1/2*ALPHA*(L-z))*ho*sin(M_PI*x/a)*sinh(beta*z)/(2*sinh(beta*L)) ;
    
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
