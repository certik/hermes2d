#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"

#include "hermes2d.h"
#include "solver_umfpack.h"

using namespace RefinementSelectors;

// This example uses automatic adaptivity to solve a 4-group neutron diffusion equation in the reactor core.
// The eigenproblem is solved using power interations.
//
// The reactor neutronics in a general coordinate system is given by the following eigenproblem:
//
//  - \nabla \cdot D_g \nabla \phi_g + \Sigma_{Rg}\phi_g - \sum_{g' \neq g} \Sigma_s^{g'\to g} \phi_{g'} =
//  = \frac{\chi_g}{k_{eff}} \sum_{g'} \nu_{g'} \Sigma_{fg'}\phi_{g'}
//
// where 1/k_{eff} is eigenvalue and \phi_g, g = 1,...,4 are eigenvectors (neutron fluxes). The current problem
// is posed in a 3D cylindrical axisymmetric geometry, leading to a 2D problem with r-z as the independent spatial 
// coordinates. The corresponding diffusion operator is given by (r = x, z = y):
//
//	\nabla \cdot D \nabla \phi = \frac{1}{x} (x D \phi_x)_x  + (D \phi_y)_y 
//
// BC:
//
// homogeneous neumann on symmetry axis
// d D_g\phi_g / d n = - 0.5 \phi_g   elsewhere
//
// The eigenproblem is numerically solved using common technique known as the power method (power iterations):
//
//  1) Make an initial estimate of \phi_g and k_{eff}
//  2) For n = 1, 2,...
//         solve for \phi_g using previous k_prev
//         solve for new k_{eff}
//                                \int_{Active Core} \sum^4_{g = 1} \nu_{g} \Sigma_{fg}\phi_{g}_{new}
//               k_new =  k_prev -------------------------------------------------------------------------
//                                \int_{Active Core} \sum^4_{g = 1} \nu_{g} \Sigma_{fg}\phi_{g}_{prev}
//  3) Stop iterations when
//
//     |   k_new - k_prev  |
//     | ----------------- |  < epsilon
//     |       k_new       |
//
//  The following parameters can be changed:

#include "h1_adapt_custom_norm.h"

const bool SOLVE_ON_COARSE_MESH = false; // If true, coarse mesh FE problem is solved in every adaptivity step.
                                         // If false, projection of the fine mesh solution on the coarse mesh is used. 
const int INIT_REF_NUM[4] = {1, 1, 1, 1};// Initial uniform mesh refinement for the individual solution components.
const int P_INIT[4] = {1, 1, 1, 1};      // Initial polynomial orders for the individual solution components.	
const double THRESHOLD = 0.3;            // This is a quantitative parameter of the adapt(...) function and
                                         // it has different meanings for various adaptive strategies (see below).
const int STRATEGY = 1;                  // Adaptive strategy:
                                         // STRATEGY = 0 ... refine elements until sqrt(THRESHOLD) times total
                                         //   error is processed. If more elements have similar errors, refine
                                         //   all to keep the mesh symmetric.
                                         // STRATEGY = 1 ... refine all elements whose error is larger
                                         //   than THRESHOLD times maximum element error.
                                         // STRATEGY = 2 ... refine all elements whose error is larger
                                         //   than THRESHOLD.
                                         // More adaptive strategies can be created in adapt_ortho_h1.cpp.
const CandList CAND_LIST = H2D_HP_ANISO; // Predefined list of element refinement candidates. Possible values are
                                         // H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO,
                                         // H2D_HP_ANISO_H, H2D_HP_ANISO_P, H2D_HP_ANISO.
                                         // See User Documentation for details.
const int MESH_REGULARITY = -1;          // Maximum allowed level of hanging nodes:
                                         // MESH_REGULARITY = -1 ... arbitrary level hangning nodes (default),
                                         // MESH_REGULARITY = 1 ... at most one-level hanging nodes,
                                         // MESH_REGULARITY = 2 ... at most two-level hanging nodes, etc.
                                         // Note that regular meshes are not supported, this is due to
                                         // their notoriously bad performance.
const double CONV_EXP = 1.0;             // Default value is 1.0. This parameter influences the selection of
                                         // candidates in hp-adaptivity. See get_optimal_refinement() for details.
const double ERR_STOP = 0.01;            // Stopping criterion for adaptivity (rel. error tolerance between the
                                         // fine mesh and coarse mesh solution in percent).
const int NDOF_STOP = 60000;             // Adaptivity process stops when the number of degrees of freedom grows over
                                         // this limit. This is mainly to prevent h-adaptivity to go on forever.
const int MAX_ADAPT_NUM = 30;	         	 // Adaptivity process stops when the number of adaptation steps grows over
                                         // this limit.
const NormType NORM_TYPE = NORM_DEFAULT; // Specifies norm used by the class H1AdaptCustomNorm to calculate relative error 
																				 // NORM_DEFAULT = 0 ... the same norm as used by H1Adapt 
																				 //		(defined either by the supplied bilinear forms or by the standard H1 norm)
																				 // NORM_MAX = 1 ... maximum (L_infty) norm


// Macro for simpler definition of bilinear forms in the energy norm.
#define callback_egnorm(a)     a<scalar, scalar>, a<Ord, Ord>

// Element markers.
const int marker_reflector = 1;
const int marker_core = 2;

// Boundary markers.
const int bc_vacuum = 1;
const int bc_sym = 2;

// Boundary condition types.
BCType bc_types(int marker)
{
  return BC_NATURAL;
}

// Essential (Dirichlet) boundary condition values.
scalar essential_bc_values(int ess_bdy_marker, double x, double y)
{
  return 0;
}

// Reflector properties (0) core properties (1),
const double D[2][4] = {{0.0164, 0.0085, 0.00832, 0.00821},
                        {0.0235, 0.0121, 0.0119, 0.0116}};
const double Sa[2][4] = {{0.00139, 0.000218, 0.00197, 0.0106},
                         {0.00977, 0.162, 0.156, 0.535}};
const double Sr[2][4] = {{1.77139, 0.533218, 3.31197, 0.0106},
                         {1.23977, 0.529, 2.436, 0.535}};
const double Sf[2][4] = {{0.0, 0.0, 0.0, 0.0}, {0.00395, 0.0262, 0.0718, 0.346}};
const double nu[2][4] = {{0.0, 0.0, 0.0, 0.0}, {2.49, 2.43, 2.42, 2.42}};
const double chi[2][4] = {{0.0, 0.0, 0.0, 0.0}, {0.9675, 0.03250, 0.0, 0.0}};
const double Ss[2][4][4] = {{{ 0.0,   0.0,  0.0, 0.0},
                             {1.77,   0.0,  0.0, 0.0},
                             { 0.0, 0.533,  0.0, 0.0},
                             { 0.0,   0.0, 3.31, 0.0}},
                            {{ 0.0,   0.0,  0.0, 0.0},
                             {1.23,   0.0,  0.0, 0.0},
                             { 0.0, 0.367,  0.0, 0.0},
                             { 0.0,   0.0, 2.28, 0.0}}};

// Power iteration control.												
double k_eff = 1.0;			// Initial eigenvalue approximation.
double TOL_PIT_CM = 5e-5;		// Tolerance for eigenvalue convergence when solving on coarse mesh.
double TOL_PIT_RM = 5e-7;		// Tolerance for eigenvalue convergence when solving on reference mesh.

// Weak forms.
#include "forms.cpp"

// Source function.
void source_fn(int n, scalar* a, scalar* b, scalar* c, scalar* d, scalar* out)
{
  for (int i = 0; i < n; i++)
  {
    out[i] = (nu[1][0] * Sf[1][0] * a[i] +
        nu[1][1] * Sf[1][1] * b[i] +
        nu[1][2] * Sf[1][2] * c[i] +
        nu[1][3] * Sf[1][3] * d[i]);
  }
}

// Integral over the active core.
double integrate(MeshFunction* sln, int marker)
{
  Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);

  double integral = 0.0;
  Element* e;
  Mesh* mesh = sln->get_mesh();

  for_all_active_elements(e, mesh)
  {
    if (e->marker == marker)
    {
      update_limit_table(e->get_mode());
      sln->set_active_element(e);
      RefMap* ru = sln->get_refmap();
      int o = sln->get_fn_order() + ru->get_inv_ref_order();
      limit_order(o);
      sln->set_quad_order(o, H2D_FN_VAL);
      scalar *uval = sln->get_fn_values();
      double* x = ru->get_phys_x(o);
      double result = 0.0;
      h1_integrate_expression(x[i] * uval[i]);
      integral += result;
    }
  }

  return 2.0 * M_PI * integral;
}

// Calculate number of negative solution values.
int get_num_of_neg(MeshFunction *sln)
{
	Quad2D* quad = &g_quad_2d_std;
  sln->set_quad_2d(quad);
  Element* e;
  Mesh* mesh = sln->get_mesh();
  
  int n = 0;

  for_all_active_elements(e, mesh)
  {
    update_limit_table(e->get_mode());
    sln->set_active_element(e);
    RefMap* ru = sln->get_refmap();
    int o = sln->get_fn_order() + ru->get_inv_ref_order();
    limit_order(o);
    scalar *uval = sln->get_fn_values();
    int np = quad->get_num_points(o);
			
		for (int i = 0; i < np; i++)
			if (uval[i] < -1e-12)
				n++;
  }
  
  return n;
}

// Power iteration using "ls" as the linear system associated with the generalized 
// eigenvalue problem and "iter_X" as the initial guess for eigenvector; "ls" 
// is assumed to be already assembled and only rhs updates are performed; converged
// eigenvectors are stored in "sln_X" and the eigenvalue in global variable "k_eff".
void power_iteration(Solution *sln1, Solution *sln2, Solution *sln3, Solution *sln4,
		     Solution *iter1, Solution *iter2, Solution *iter3, Solution *iter4,
		     LinSystem *ls, double tol)
{
  bool eigen_done = false; int it = 0;
  do {
    // Solve for new eigenvectors.
    ls->solve(Tuple<Solution*>(sln1, sln2, sln3, sln4));

    // Update fission sources.
    SimpleFilter source(source_fn, sln1, sln2, sln3, sln4);
    SimpleFilter source_prev(source_fn, iter1, iter2, iter3, iter4);

    // Compute eigenvalue.
    double k_new = k_eff * (integrate(&source, marker_core) / integrate(&source_prev, marker_core));
	  
    info("      dominant eigenvalue (est): %g, rel error: %g", k_new, fabs((k_eff - k_new) / k_new));

    // Stopping criterion.
    if (fabs((k_eff - k_new) / k_new) < tol) eigen_done = true;

    // Store eigenvectors for next iteration.
    iter1->copy(sln1);    
    iter2->copy(sln2);
    iter3->copy(sln3);    
    iter4->copy(sln4);
	  
    // Update eigenvalue.
    k_eff = k_new;
	  	  
	if (!eigen_done)
		// Update rhs of the system considering the updated eigenpair approximation.
		ls->assemble(true);
	}
	while (!eigen_done);
}

int main(int argc, char* argv[])
{
  // Time measurement.
  TimePeriod cpu_time;
  cpu_time.tick();

  // Load the mesh.
  Mesh mesh1, mesh2, mesh3, mesh4;
  H2DReader mloader;
  mloader.load("reactor.mesh", &mesh1);

  // Obtain meshes for the 2nd to 4th group by cloning the mesh loaded for the 1st group.
  // This initializes the multimesh hp-FEM.
  mesh2.copy(&mesh1);
  mesh3.copy(&mesh1);
  mesh4.copy(&mesh1);

  // Initial uniform refinements.
  info("Setting initial conditions.");
  for (int i = 0; i < INIT_REF_NUM[0]; i++) mesh1.refine_all_elements();
  for (int i = 0; i < INIT_REF_NUM[1]; i++) mesh2.refine_all_elements();
  for (int i = 0; i < INIT_REF_NUM[2]; i++) mesh3.refine_all_elements();
  for (int i = 0; i < INIT_REF_NUM[3]; i++) mesh4.refine_all_elements();

  // Solution variables.
  Solution iter1, iter2, iter3, iter4,                          // Previous iterations.
           sln1_coarse, sln2_coarse, sln3_coarse, sln4_coarse,	// Current iterations on coarse meshes.
	   sln1_fine, sln2_fine, sln3_fine, sln4_fine;	        // Current iterations on fine meshes.

  // Set initial conditions for the Newton's method.
  iter1.set_const(&mesh1, 1.00);
  iter2.set_const(&mesh2, 1.00);
  iter3.set_const(&mesh3, 1.00);
  iter4.set_const(&mesh4, 1.00);

  // Matrix solver.
  UmfpackSolver umfpack;

  // Create H1 spaces with default shapesets.
  H1Space space1(&mesh1, bc_types, essential_bc_values, P_INIT[0]);
  H1Space space2(&mesh2, bc_types, essential_bc_values, P_INIT[1]); 
  H1Space space3(&mesh3, bc_types, essential_bc_values, P_INIT[2]); 
  H1Space space4(&mesh4, bc_types, essential_bc_values, P_INIT[3]); 

  // Initialize the weak formulation.
  WeakForm wf(4);
  wf.add_matrix_form(0, 0, callback(biform_0_0), H2D_SYM);
  wf.add_matrix_form(1, 1, callback(biform_1_1), H2D_SYM);
  wf.add_matrix_form(1, 0, callback(biform_1_0));
  wf.add_matrix_form(2, 2, callback(biform_2_2), H2D_SYM);
  wf.add_matrix_form(2, 1, callback(biform_2_1));
  wf.add_matrix_form(3, 3, callback(biform_3_3), H2D_SYM);
  wf.add_matrix_form(3, 2, callback(biform_3_2));
  wf.add_vector_form(0, callback(liform_0), marker_core, Tuple<MeshFunction*>(&iter1, &iter2, &iter3, &iter4));
  wf.add_vector_form(1, callback(liform_1), marker_core, Tuple<MeshFunction*>(&iter1, &iter2, &iter3, &iter4));
  wf.add_vector_form(2, callback(liform_2), marker_core, Tuple<MeshFunction*>(&iter1, &iter2, &iter3, &iter4));
  wf.add_vector_form(3, callback(liform_3), marker_core, Tuple<MeshFunction*>(&iter1, &iter2, &iter3, &iter4));
  wf.add_matrix_form_surf(0, 0, callback(biform_surf_0_0), bc_vacuum);
  wf.add_matrix_form_surf(1, 1, callback(biform_surf_1_1), bc_vacuum);
  wf.add_matrix_form_surf(2, 2, callback(biform_surf_2_2), bc_vacuum);
  wf.add_matrix_form_surf(3, 3, callback(biform_surf_3_3), bc_vacuum);

  // Initialize and solve coarse mesh problem.
  LinSystem ls(&wf, &umfpack, Tuple<Space*>(&space1, &space2, &space3, &space4));
  ls.assemble();
  info("Coarse mesh power iteration, %d + %d + %d + %d = %d ndof:", 
  ls.get_num_dofs(0), ls.get_num_dofs(1), ls.get_num_dofs(2), ls.get_num_dofs(3), ls.get_num_dofs()); 
  power_iteration(&sln1_coarse, &sln2_coarse, &sln3_coarse, &sln4_coarse, 
	          &iter1, &iter2, &iter3, &iter4,
		  &ls, TOL_PIT_CM);

  // Initialize views.
  ScalarView view1("Neutron flux 1", 0, 0, 320, 400);
  ScalarView view2("Neutron flux 2", 330, 0, 320, 400);
  ScalarView view3("Neutron flux 3", 660, 0, 320, 400);
  ScalarView view4("Neutron flux 4", 990, 0, 320, 400);
  OrderView oview1("Mesh for group 1", 0, 450, 320, 500);
  OrderView oview2("Mesh for group 2", 330, 450, 320, 500);
  OrderView oview3("Mesh for group 3", 660, 450, 320, 500);
  OrderView oview4("Mesh for group 4", 990, 450, 320, 500);

  // Do not show meshes.
  view1.show_mesh(false);
  view2.show_mesh(false);
  view3.show_mesh(false);
  view4.show_mesh(false);
  
  // DOF and CPU convergence graphs.
  SimpleGraph graph_dof, graph_cpu, graph_dof_keff, graph_cpu_keff;

  // Initialize refinement selector.
  H1ProjBasedSelector selector(CAND_LIST, CONV_EXP, H2DRS_DEFAULT_ORDER);

  // Adaptivity loop:
  int as = 1; bool done = false;
  do {
		
    info("---- Adaptivity step %d:", as);

    // Initialize reference problem. 
    RefSystem rs(&ls);

    // First time project coarse mesh solutions on fine meshes.
    if (as == 1) {
      info("Projecting coarse mesh solutions on fine meshes.");
      rs.project_global(Tuple<MeshFunction*>(&sln1_coarse, &sln2_coarse, &sln3_coarse, &sln4_coarse), 
                        Tuple<Solution*>(&iter1, &iter2, &iter3, &iter4));
    }

    // Solve the fine mesh problem.
    rs.assemble();	
    info("Fine mesh power iteration, %d + %d + %d + %d = %d ndof:", 
      rs.get_num_dofs(0), rs.get_num_dofs(1), rs.get_num_dofs(2), rs.get_num_dofs(3), rs.get_num_dofs());
    power_iteration(&sln1_fine, &sln2_fine, &sln3_fine, &sln4_fine,
		    &iter1, &iter2, &iter3, &iter4,
		    &rs, TOL_PIT_RM);

    // Either solve on coarse mesh or project the fine mesh solution 
    // on the coarse mesh.
    if (SOLVE_ON_COARSE_MESH) {
      if (as > 1) {
        info("Coarse mesh power iteration, %d + %d + %d + %d = %d ndof:", 
        ls.get_num_dofs(0), ls.get_num_dofs(1), ls.get_num_dofs(2), ls.get_num_dofs(3), ls.get_num_dofs()); 
        ls.assemble();
        power_iteration(&sln1_coarse, &sln2_coarse, &sln3_coarse, &sln4_coarse, 
	  	        &iter1, &iter2, &iter3, &iter4,
		        &ls, TOL_PIT_CM);
      }
    }
    else {
      info("Projecting fine mesh solutions on coarse meshes.");
      ls.project_global(Tuple<MeshFunction*>(&sln1_fine, &sln2_fine, &sln3_fine, &sln4_fine), 
                        Tuple<Solution*>(&sln1_coarse, &sln2_coarse, &sln3_coarse, &sln4_coarse));
    }

    // Time measurement.
    cpu_time.tick();

    // View the solution and meshes.
    view1.show(&sln1_coarse);
    view2.show(&sln2_coarse);
    view3.show(&sln3_coarse);
    view4.show(&sln4_coarse);
    oview1.show(&space1);
    oview2.show(&space2);
    oview3.show(&space3);
    oview4.show(&space4);

    // Skip visualization time.
    cpu_time.tick(H2D_SKIP);
    
    // Report the number of negative eigenfunction values.
    info("Num. of negative values: %d, %d, %d, %d", 
	 get_num_of_neg(&sln1_coarse), get_num_of_neg(&sln2_coarse), 
         get_num_of_neg(&sln3_coarse), get_num_of_neg(&sln4_coarse));		
		    
    // Calculate element errors and total error estimate.
    H1AdaptCustomNorm hp(NORM_TYPE, Tuple<Space*>(&space1, &space2, &space3, &space4));
    hp.set_biform(0, 0, callback_egnorm(biform_0_0));
  	hp.set_biform(1, 1, callback_egnorm(biform_1_1));
		hp.set_biform(1, 0, callback_egnorm(biform_1_0));
		hp.set_biform(2, 2, callback_egnorm(biform_2_2));
		hp.set_biform(2, 1, callback_egnorm(biform_2_1));
		hp.set_biform(3, 3, callback_egnorm(biform_3_3));
		hp.set_biform(3, 2, callback_egnorm(biform_3_2));
		
    // Calculate element errors and error estimate for adaptivity.
    info("Calculating error.");
    hp.set_solutions(Tuple<Solution*>(&sln1_coarse, &sln2_coarse, &sln3_coarse, &sln4_coarse),
                     Tuple<Solution*>(&sln1_fine, &sln2_fine, &sln3_fine, &sln4_fine));
    double err_est = hp.calc_error(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_REL) * 100;
  
    // Time measurement.
    cpu_time.tick();        
    double cta = cpu_time.accumulated();

    // Calculate H1 error estimates.
    double err_est_1 = h1_error(&sln1_coarse, &sln1_fine);
    double err_est_2 = h1_error(&sln2_coarse, &sln2_fine);
    double err_est_3 = h1_error(&sln3_coarse, &sln3_fine);
    double err_est_4 = h1_error(&sln4_coarse, &sln4_fine);

    // Report results.
    info("ndof_coarse: %d + %d + %d + %d = %d", ls.get_num_dofs(0), ls.get_num_dofs(1), ls.get_num_dofs(2), 
	 ls.get_num_dofs(3), ls.get_num_dofs());  
    info("err_est_coarse: %g%%, %g%%, %g%%, %g%%", err_est_1, err_est_2, err_est_3, err_est_4); 
  
    // eigenvalue error w.r.t. solution obtained on a 3x uniformly refined mesh
  	// with uniform distribution of polynomial degrees (=4), converged to within
  	// tolerance of 5e-11; in units of percent-milli (pcm)
  	double keff_err = 1e5*fabs(k_eff - 1.14091443)/1.14091443;
  					
    // Add entry to DOF convergence graph.
    graph_dof.add_values(ls.get_num_dofs(), err_est);
    graph_dof.save("conv_dof.dat");
  
    // Add entry to CPU convergence graph.
    graph_cpu.add_values(cta, err_est);
    graph_cpu.save("conv_cpu.dat");
    
    // Add entry to DOF convergence graph w.r.t. dominant eigenvalue.
    graph_dof_keff.add_values(ls.get_num_dofs(), keff_err);
    graph_dof_keff.save("conv_dof_keff.dat");
    
    // Add entry to CPU convergence graph w.r.t. dominant eigenvalue.
    graph_cpu_keff.add_values(cta, keff_err);
    graph_cpu_keff.save("conv_cpu_keff.dat");
    
    // If err_est too large, adapt the mesh.
    if (err_est < ERR_STOP) break;
    else {
      info("Adapting the coarse mesh.");
      done = hp.adapt(&selector, THRESHOLD, STRATEGY, MESH_REGULARITY);
      if (ls.get_num_dofs() >= NDOF_STOP) done = true;
    }

    as++;
    if (as >= MAX_ADAPT_NUM) done = true;
  }
  while(done == false);	
  verbose("Total running time: %g s", cpu_time.accumulated());
  
  // Wait for all views to be closed.
  View::wait();
  return 0;
};
