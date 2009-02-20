void calculate_elements_length(double* ele_len, double* u_infty, Solution* u1, Solution* u2, Mesh* mesh)
{
  Quad2D* quad = &g_quad_2d_std;
  u1->set_quad_2d(quad);
  u2->set_quad_2d(quad);
  
  Solution tmp;
  tmp.set_zero(mesh);

  Mesh* meshes[3] = { mesh, u1->get_mesh(), u2->get_mesh() };
  Transformable* tr[3] = { &tmp, u1, u2 };
  Traverse trav;
  trav.begin(3, meshes, tr);
  
  int ne = mesh->get_max_element_id() + 1;
  double* a = new double[ne];
  double* b = new double[ne];
  int* n = new int[ne];
  memset(a, 0, sizeof(double) * ne);
  memset(b, 0, sizeof(double) * ne);
  memset(n, 0, sizeof(int) * ne);

  Element** ee; 
  while ((ee = trav.get_next_state(NULL, NULL)) != NULL)
  {
    int o = u1->get_fn_order();
    u1->set_quad_order(o);
    u2->set_quad_order(o);
    scalar *uval, *vval;
    uval = u1->get_fn_values();
    vval = u2->get_fn_values();
    int np = quad->get_num_points(o);
    for (int i = 0; i < np; i++) {
      a[ee[0]->id] += uval[i];
      b[ee[0]->id] += vval[i];
    }
    n[ee[0]->id] += np;
  }
  trav.finish();

  Element* e;
  for_all_active_elements(e, mesh)
  {
    // averaged values of velocities on elements of mesh
    double c = a[e->id] / n[e->id]; 
    double d = b[e->id] / n[e->id]; 
    u_infty[e->id] = (fabs(c) > fabs(d)) ? fabs(c) : fabs(d);

    double den = (sqrt(sqr(c) + sqr(d)));
    // find the max length
    double min = 10000.0, max = -10000.0;
    for (int i = 0; i < e->nvert; i++)
    {
      double x = e->vn[i]->x;
      double y = e->vn[i]->y;
      double l = (c*x + d*y) / den;
      if (l > max) max = l;
      if (l < min) min = l;
    }

    ele_len[e->id] = max - min;
  }

  delete [] a;
  delete [] b;
  delete [] n;

}


void calculate_stabilization_parameters(double* delta_K, double* tau_K, Solution* xsln, Solution* ysln, Mesh* mesh)
{
  int ne = mesh->get_max_element_id() + 1;
  double* ele_len = new double[ne];
  double* u_infty = new double[ne];
  
  calculate_elements_length(ele_len, u_infty, xsln, ysln, mesh); 

  Element* e;
  for_all_active_elements(e, mesh)
  {
    delta_K[e->id] = delta_star * sqr(ele_len[e->id]);
    tau_K[e->id] = tau_star * 1.0;
  }

  delete [] ele_len;
  delete [] u_infty;
}



