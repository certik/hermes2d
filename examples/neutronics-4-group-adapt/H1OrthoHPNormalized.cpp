// function used to calculate error in H1 norm
double error_fn_h1(MeshFunction* sln1, MeshFunction* sln2, RefMap* ru, RefMap* rv)
{
  Quad2D* quad = sln1->get_quad_2d();

  int o = 2*std::max(sln1->get_fn_order(), sln2->get_fn_order()) + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln1->set_quad_order(o);
  sln2->set_quad_order(o);

  scalar *uval, *vval, *dudx, *dudy, *dvdx, *dvdy;
  uval = sln1->get_fn_values();
  vval = sln2->get_fn_values();
  sln1->get_dx_dy_values(dudx, dudy);
  sln2->get_dx_dy_values(dvdx, dvdy);

  double result = 0.0;
  h1_integrate_expression(sqr(uval[i] - vval[i]) +
                          sqr(dudx[i] - dvdx[i]) + sqr(dudy[i] - dvdy[i]));
  return result;
}

// function used to calculate H1 norm of the solution
double norm_fn_h1(MeshFunction* sln, RefMap* ru)
{
  Quad2D* quad = sln->get_quad_2d();

  int o = 2 * sln->get_fn_order() + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln->set_quad_order(o);

  scalar *uval, *dudx, *dudy;
  uval = sln->get_fn_values();
  sln->get_dx_dy_values(dudx, dudy);

  double result = 0.0;
  h1_integrate_expression(sqr(uval[i]) + sqr(dudx[i]) + sqr(dudy[i]));
  return result;
}

// function used to calculate L_oo norm of the solution
double norm_fn_inf(MeshFunction* sln, RefMap* ru)
{
  Quad2D* quad = sln->get_quad_2d();

  int o = 2 * sln->get_fn_order() + ru->get_inv_ref_order();
  limit_order_nowarn(o);

  sln->set_quad_order(o);

  scalar *uval;
  uval = sln->get_fn_values();

  double result = uval[0];
  int np = quad->get_num_points(o);
  for (int i = 0; i < np; i++)
	  if (uval[i] > result)
	  	result = uval[i];
	  	
  return result;
}

class H2D_API H1OrthoHPNormalized : public H1OrthoHP
{
	public:
		H1OrthoHPNormalized(int num, Space* space1, Space* space2 = NULL, Space* space3 = NULL, Space* space4 = NULL, Space* space5 = NULL, Space* space6 = NULL, Space* space7 = NULL, Space* space8 = NULL, Space* space9 = NULL, Space* space10 = NULL)
		: H1OrthoHP(num, space1, space2, space3, space4, space5, space6, space7, space8, space9, space10) {};
		
		void sort_elements_by_error(Mesh** meshes);
		
	private:
  	static double** cmp_err; ///< An helper array used to sort reference to elements.
  	static int compare(const void* p1, const void* p2); ///< Compares to reference to an element according to H1OrthoHP::cmp_err.

};

class H2D_API H1OrthoHPNormalized2 : public H1OrthoHPNormalized
{
	public:
		H1OrthoHPNormalized2(int norm, int num, Space* space1, Space* space2 = NULL, Space* space3 = NULL, Space* space4 = NULL, Space* space5 = NULL, Space* space6 = NULL, Space* space7 = NULL, Space* space8 = NULL, Space* space9 = NULL, Space* space10 = NULL)
		: H1OrthoHPNormalized(num, space1, space2, space3, space4, space5, space6, space7, space8, space9, space10) { norm_type = norm; }
		
		double calc_error_n(int n, ...);
		
	private:
		int norm_type;
};

double** H1OrthoHPNormalized::cmp_err;
int H1OrthoHPNormalized::compare(const void* p1, const void* p2) {
  const ElementReference& e1 = *((const ElementReference*)p1);
  const ElementReference& e2 = *((const ElementReference*)p2);
  return cmp_err[e1.comp][e1.id] < cmp_err[e2.comp][e2.id] ? 1 : -1;
}

void H1OrthoHPNormalized::sort_elements_by_error(Mesh** meshes) {
  //allocate
  if (esort != NULL)
    delete[] esort;
  esort = new ElementReference[nact];

  //prepare indices
  Element* e;
  int inx = 0;
  for (int i = 0; i < num; i++)
    for_all_active_elements(e, meshes[i]) {
      esort[inx].id = e->id;
      esort[inx].comp = i;
      inx++;
      errors[i][e->id] = errors[i][e->id] / norms[i]; // sqrt or not?
    }

  //sort
  assert(inx == nact);
  cmp_err = errors;
  qsort(esort, nact, sizeof(ElementReference), compare);
}

double H1OrthoHPNormalized2::calc_error_n(int n, ...)
{
  int i, j;

  if (n != num) error("Wrong number of solutions.");
  
  // obtain solutions and bilinear forms
  va_list ap;
  va_start(ap, n);
  for (i = 0; i < n; i++) {
    sln[i] = va_arg(ap, Solution*); //?WTF: input of calc_error, which calls calc_error_n, is a type MeshFunction* that is parent of Solution*
    sln[i]->set_quad_2d(&g_quad_2d_std);
  }
  for (i = 0; i < n; i++) {
    rsln[i] = va_arg(ap, Solution*); //?WTF: input of calc_error, which calls calc_error_n, is a type MeshFunction* that is parent of Solution*
    rsln[i]->set_quad_2d(&g_quad_2d_std);
  }
  va_end(ap);
  
  // prepare multi-mesh traversal and error arrays
  AUTOLA_OR(Mesh*, meshes, 2*num);
  AUTOLA_OR(Transformable*, tr, 2*num);
  Traverse trav;
	nact = 0;
	
  for (i = 0; i < num; i++)
  {
    meshes[i] = sln[i]->get_mesh();
    meshes[i+num] = rsln[i]->get_mesh();
    tr[i] = sln[i];
    tr[i+num] = rsln[i];

		nact += sln[i]->get_mesh()->get_num_active_elements();

    int max = meshes[i]->get_max_element_id();
    if (errors[i] != NULL) delete [] errors[i];
    errors[i] = new double[max];
    memset(errors[i], 0, sizeof(double) * max);
  }

  double total_norm = 0.0;
  AUTOLA_OR(double, norms, num);
  memset(norms, 0, norms.size);
  double total_error = 0.0;
  
  Element** ee;
  trav.begin(2*num, meshes, tr);
  while ((ee = trav.get_next_state(NULL, NULL)) != NULL)
  {
    for (i = 0; i < num; i++)
    {
      RefMap* rm = sln[i]->get_refmap();
      RefMap* rrm = rsln[i]->get_refmap();

      double e, t;
      e = error_fn_h1(sln[i], rsln[i], rm, rrm);
      
      if (norm_type == 0)
	      t = norm_fn_h1(rsln[i], rrm);
	    else
	    	t = norm_fn_inf(rsln[i], rrm);
      
      norms[i] += t;
      total_norm += t;
      total_error += e;
      errors[i][ee[i]->id] += e;
    }
  }
  trav.finish();  
   
  //prepare an ordered list of elements according to an error
  sort_elements_by_error(meshes);

  have_errors = true;
  total_err = total_error / total_norm;
  return sqrt(total_error / total_norm);
}
