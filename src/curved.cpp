// This file is part of Hermes2D.
//
// Hermes2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Hermes2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Hermes2D.  If not, see <http://www.gnu.org/licenses/>.

#include "common.h"
#include "shapeset_h1_all.h"
#include "shapeset_common.h"
#include "precalc.h"
#include "curved.h"
#include "mesh.h"
#include "quad_all.h"
#include "matrix_old.h"


// defined in refmap.cpp
extern H1ShapesetBeuchler ref_map_shapeset;
extern PrecalcShapeset ref_map_pss;

static double** edge_proj_matrix = NULL;  //projection matrix for each edge is the same
static double** bubble_proj_matrix_tri = NULL; //projection matrix for triangle bubbles
static double** bubble_proj_matrix_quad = NULL; //projection matrix for quad bubbles

static double* edge_p = NULL;  // diagonal vector in cholesky factorization
static double* bubble_tri_p = NULL; // diagonal vector in cholesky factorization
static double* bubble_quad_p = NULL; // diagonal vector in cholesky factorization

static Quad1DStd quad1d;
static Quad2DStd quad2d; // fixme: g_quad_2d_std

static Trf ctm;


//// NURBS //////////////////////////////////////////////////////////////////////////////////////////

// recursive calculation of the basis function N_i,k
double nurbs_basis_fn(int i, int k, double t, double* knot)
{
  if (k == 0)
  {
    return (t >= knot[i] && t <= knot[i+1] && knot[i] < knot[i+1]) ? 1.0 : 0.0;
  }
  else
  {
    double N1 = nurbs_basis_fn(i, k-1, t, knot);
    double N2 = nurbs_basis_fn(i+1, k-1, t, knot);

    double result = 0.0;
    if (knot[i+k] != knot[i])
    {
      result += ((t - knot[i]) / (knot[i+k] - knot[i])) * N1;
    }
    if (knot[i+k+1] != knot[i+1])
    {
      result += ((knot[i+k+1] - t) / (knot[i+k+1] - knot[i+1])) * N2;
    }
    return result;
  }
}


// nurbs curve: t goes from -1 to 1, function returns x, y coordinates in plane
void nurbs_edge(Element* e, Nurbs* nurbs, int edge, double t, double& x, double& y)
{
  t = (t + 1) / 2.0; // nurbs curves are parametrized from 0 to 1
  if (nurbs == NULL)
  {
    double2 v;
    v[0] = e->vn[e->next_vert(edge)]->x - e->vn[edge]->x;
    v[1] = e->vn[e->next_vert(edge)]->y - e->vn[edge]->y;
    x = e->vn[edge]->x + t * v[0];
    y = e->vn[edge]->y + t * v[1];
  }
  else
  {
    double3* cp = nurbs->pt;
    x = y = 0.0;
    double sum = 0.0;  // sum of basis fns and weights

    for (int i = 0; i < nurbs->np; i++)
    {
      double basis = nurbs_basis_fn(i, nurbs->degree, t, nurbs->kv);
      sum += cp[i][2] * basis;
      x   += cp[i][2] * basis * cp[i][0];
      y   += cp[i][2] * basis * cp[i][1];
    }

    sum = 1.0 / sum;
    x *= sum;
    y *= sum;
  }
}



//// non-polynomial reference map //////////////////////////////////////////////////////////////////////////////////

// definition of vertex basis functions for triangle
static double lambda_0(double x, double y) { return -0.5 * (x + y); }
static double lambda_1(double x, double y) { return  0.5 * (x + 1); }
static double lambda_2(double x, double y) { return  0.5 * (y + 1); }

static double (*lambda[3])(double, double) = { lambda_0, lambda_1, lambda_2 };

// 1D Lobatto functions
static double lob0(double x)  { return l0(x); }
static double lob1(double x)  { return l1(x); }
static double lob2(double x)  { return l2(x); }
static double lob3(double x)  { return l3(x); }
static double lob4(double x)  { return l4(x); }
static double lob5(double x)  { return l5(x); }
static double lob6(double x)  { return l6(x); }
static double lob7(double x)  { return l7(x); }
static double lob8(double x)  { return l8(x); }
static double lob9(double x)  { return l9(x); }
static double lob10(double x) { return l10(x); }
static double lob11(double x) { return l11(x); }

static double (*lob[12])(double) = { lob0, lob1, lob2, lob3, lob4, lob5, lob6, lob7, lob8, lob9, lob10, lob11 };

static double2 ref_vert[2][4] =
{
  { { -1.0, -1.0 }, { 1.0, -1.0 }, { -1.0, 1.0 }, {  0.0, 0.0 } },
  { { -1.0, -1.0 }, { 1.0, -1.0 }, {  1.0, 1.0 }, { -1.0, 1.0 } }
};


// subtraction of straight edge and nurbs curve
static void nurbs_edge_0(Element* e, Nurbs* nurbs, int edge, double t, double& x, double& y)
{
  int va = edge;
  int vb = e->next_vert(edge);
  nurbs_edge(e, nurbs, edge, t, x, y);

  x -= 0.5 * ((1-t) * (e->vn[va]->x) + (1+t) * (e->vn[vb]->x));
  y -= 0.5 * ((1-t) * (e->vn[va]->y) + (1+t) * (e->vn[vb]->y));

  double k = 4.0 / ((1-t) * (1+t));
  x *= k;
  y *= k;
}


// calculation of nonpolynomial reference mapping on curved element
static void calc_ref_map_tri(Element* e, Nurbs** nurbs, double xi_1, double xi_2, double& x, double& y)
{
  double  fx,  fy;
  x = y = 0.0;

  for (unsigned int j = 0; j < e->nvert; j++)
  {
    int va = j;
    int vb = e->next_vert(j);
    double l_a = lambda[va](xi_1, xi_2);
    double l_b = lambda[vb](xi_1, xi_2);

    // vertex part
    x += e->vn[j]->x * l_a;
    y += e->vn[j]->y * l_a;

    if (!(((ref_vert[0][va][0] == xi_1) && (ref_vert[0][va][1] == xi_2)) ||
          ((ref_vert[0][vb][0] == xi_1) && (ref_vert[0][vb][1] == xi_2))))
    {
      // edge part
      double t = l_b - l_a;
      nurbs_edge_0(e, nurbs[j], j, t, fx, fy);
      x += fx * l_a  * l_b;
      y += fy * l_a  * l_b;
    }
  }
}


static void calc_ref_map_quad(Element* e, Nurbs** nurbs, double xi_1, double xi_2,
                              double& x, double& y)
{
  double ex[4], ey[4];

  nurbs_edge(e, nurbs[0], 0,  xi_1, ex[0], ey[0]);
  nurbs_edge(e, nurbs[1], 1,  xi_2, ex[1], ey[1]);
  nurbs_edge(e, nurbs[2], 2, -xi_1, ex[2], ey[2]);
  nurbs_edge(e, nurbs[3], 3, -xi_2, ex[3], ey[3]);

  x = (1-xi_2)/2.0 * ex[0] + (1+xi_1)/2.0 * ex[1] +
      (1+xi_2)/2.0 * ex[2] + (1-xi_1)/2.0 * ex[3] -
      (1-xi_1)*(1-xi_2)/4.0 * e->vn[0]->x - (1+xi_1)*(1-xi_2)/4.0 * e->vn[1]->x -
      (1+xi_1)*(1+xi_2)/4.0 * e->vn[2]->x - (1-xi_1)*(1+xi_2)/4.0 * e->vn[3]->x;

  y = (1-xi_2)/2.0 * ey[0] + (1+xi_1)/2.0 * ey[1] +
      (1+xi_2)/2.0 * ey[2] + (1-xi_1)/2.0 * ey[3] -
      (1-xi_1)*(1-xi_2)/4.0 * e->vn[0]->y - (1+xi_1)*(1-xi_2)/4.0 * e->vn[1]->y -
      (1+xi_1)*(1+xi_2)/4.0 * e->vn[2]->y - (1-xi_1)*(1+xi_2)/4.0 * e->vn[3]->y;
}


static void calc_ref_map(Element* e, Nurbs** nurbs, double xi_1, double xi_2, double2& f)
{
  if (e->get_mode() == H2D_MODE_QUAD)
    calc_ref_map_quad(e, nurbs, xi_1, xi_2, f[0], f[1]);
  else
    calc_ref_map_tri(e, nurbs, xi_1, xi_2, f[0], f[1]);
}


//// projection based interpolation ////////////////////////////////////////////////////////////////

// preparation of projection matrices, Cholesky factorization
static void precalculate_cholesky_projection_matrix_edge()
{
  int order = ref_map_shapeset.get_max_order();
  int n = order - 1; // number of edge basis functions
  edge_proj_matrix = new_matrix<double>(n, n);

  // calculate projection matrix of maximum order
  for (int i = 0; i < n; i++)
  {
    for (int j = i; j < n; j++)
    {
      int o = i + j + 4;
      double2* pt = quad1d.get_points(o);
      double val = 0.0;
      for (int k = 0; k < quad1d.get_num_points(o); k++)
      {
        double x = pt[k][0];
        double fi = lob[i+2](x);
        double fj = lob[j+2](x);
        val += pt[k][1] * (fi * fj);
      }
      edge_proj_matrix[i][j] = edge_proj_matrix[j][i] = val;
    }
  }

  // Cholesky factorization of the matrix
  edge_p = new double[n];
  choldc(edge_proj_matrix, n, edge_p);
}


// calculate the H1 seminorm products (\phi_i, \phi_j) for all 0 <= i,j < n, n is the number of bubble functions
static double** calculate_bubble_projection_matrix(int nb, int* indices)
{
  double** mat = new_matrix<double>(nb, nb);

  for (int i = 0; i < nb; i++)
  {
    for (int j = i; j < nb; j++)
    {
      int ii = indices[i], ij = indices[j];
      int o = ref_map_shapeset.get_order(ii) + ref_map_shapeset.get_order(ij);
      o = std::max(H2D_GET_V_ORDER(o), H2D_GET_H_ORDER(o));

      ref_map_pss.set_active_shape(ii);
      ref_map_pss.set_quad_order(o);
      double* fni = ref_map_pss.get_fn_values();

      ref_map_pss.set_active_shape(ij);
      ref_map_pss.set_quad_order(o);
      double* fnj = ref_map_pss.get_fn_values();

      double3* pt = quad2d.get_points(o);
      double val = 0.0;
      for (int k = 0; k < quad2d.get_num_points(o); k++)
        val += pt[k][2] * (fni[k] * fnj[k]);

      mat[i][j] = mat[j][i] = val;
    }
  }

  return mat;
}


static void precalculate_cholesky_projection_matrices_bubble()
{
  // *** triangles ***
  ref_map_pss.set_mode(H2D_MODE_TRIANGLE);
  int order = ref_map_shapeset.get_max_order();

  // calculate projection matrix of maximum order
  int nb = ref_map_shapeset.get_num_bubbles(order);
  int* indices = ref_map_shapeset.get_bubble_indices(order);
  bubble_proj_matrix_tri = calculate_bubble_projection_matrix(nb, indices);

  // cholesky factorization of the matrix
  bubble_tri_p = new double[nb];
  choldc(bubble_proj_matrix_tri, nb, bubble_tri_p);

  // *** quads ***
  ref_map_pss.set_mode(H2D_MODE_QUAD);
  order = ref_map_shapeset.get_max_order();
  order = H2D_MAKE_QUAD_ORDER(order, order);

  // calculate projection matrix of maximum order
  nb = ref_map_shapeset.get_num_bubbles(order);
  indices = ref_map_shapeset.get_bubble_indices(order);
  bubble_proj_matrix_quad = calculate_bubble_projection_matrix(nb, indices);

  // cholesky factorization of the matrix
  bubble_quad_p = new double[nb];
  choldc(bubble_proj_matrix_quad, nb, bubble_quad_p);
}


//// edge part of projection based interpolation ///////////////////////////////////////////////////

// compute point (x,y) in reference element, edge vector (v1, v2)
static void edge_coord(Element* e, int edge, double t, double2& x, double2& v)
{
  int mode = e->get_mode();
  double2 a, b;
  a[0] = ctm.m[0] * ref_vert[mode][edge][0] + ctm.t[0];
  a[1] = ctm.m[1] * ref_vert[mode][edge][1] + ctm.t[1];
  b[0] = ctm.m[0] * ref_vert[mode][e->next_vert(edge)][0] + ctm.t[0];
  b[1] = ctm.m[1] * ref_vert[mode][e->next_vert(edge)][1] + ctm.t[1];

  for (int i = 0; i < 2; i++)
  {
    v[i] = b[i] - a[i];
    x[i] = a[i] + (t+1.0)/2.0 * v[i];
  }
  double lenght = sqrt(v[0] * v[0] + v[1] * v[1]);
  v[0] /= lenght; v[1] /= lenght;
}


static void calc_edge_projection(Element* e, int edge, Nurbs** nurbs, int order, double2* proj)
{
  ref_map_pss.set_active_element(e);

  int i, j, k;
  int mo1 = quad1d.get_max_order();
  int np = quad1d.get_num_points(mo1);
  int ne = order - 1;
  int mode = e->get_mode();

  assert(np <= 15 && ne <= 10);
  double2 fn[15];
  double rhside[2][10];
  memset(fn, 0, sizeof(double2) * np);
  memset(rhside[0], 0, sizeof(double) * ne);
  memset(rhside[1], 0, sizeof(double) * ne);

  double a_1, a_2, b_1, b_2;
  a_1 = ctm.m[0] * ref_vert[mode][edge][0] + ctm.t[0];
  a_2 = ctm.m[1] * ref_vert[mode][edge][1] + ctm.t[1];
  b_1 = ctm.m[0] * ref_vert[mode][e->next_vert(edge)][0] + ctm.t[0];
  b_2 = ctm.m[1] * ref_vert[mode][e->next_vert(edge)][1] + ctm.t[1];

  // values of nonpolynomial function in two vertices
  double2 fa, fb;
  calc_ref_map(e, nurbs, a_1, a_2, fa);
  calc_ref_map(e, nurbs, b_1, b_2, fb);

  double2* pt = quad1d.get_points(mo1);
  for (j = 0; j < np; j++) // over all integration points
  {
    double2 x, v;
    double t = pt[j][0];
    edge_coord(e, edge, t, x, v);
    calc_ref_map(e, nurbs, x[0], x[1], fn[j]);

    for (k = 0; k < 2; k++)
      fn[j][k] = fn[j][k] - (fa[k] + (t+1)/2.0 * (fb[k] - fa[k]));
  }

  double2* result = proj + e->nvert + edge * (order - 1);
  for (k = 0; k < 2; k++)
  {
    for (i = 0; i < ne; i++)
    {
      for (j = 0; j < np; j++)
      {
        double t = pt[j][0];
        double fi = lob[i+2](t);
        rhside[k][i] += pt[j][1] * (fi * fn[j][k]);
      }
    }
    // solve
    cholsl(edge_proj_matrix, ne, edge_p, rhside[k], rhside[k]);
    for (i = 0; i < ne; i++)
      result[i][k] = rhside[k][i];
  }
}


//// bubble part of projection based interpolation /////////////////////////////////////////////////

static void old_projection(Element* e, int order, double2* proj, double* old[2])
{
  int mo2 = quad2d.get_max_order();
  int np = quad2d.get_num_points(mo2);

  for (unsigned int k = 0; k < e->nvert; k++) // loop over vertices
  {
    // vertex basis functions in all integration points
    double* vd;
    int index_v = ref_map_shapeset.get_vertex_index(k);
    ref_map_pss.set_active_shape(index_v);
    ref_map_pss.set_quad_order(mo2);
    vd = ref_map_pss.get_fn_values();

    for (int m = 0; m < 2; m++)   // part 0 or 1
      for (int j = 0; j < np; j++)
        old[m][j] += proj[k][m] * vd[j];

    for (int ii = 0; ii < order - 1; ii++)
    {
      // edge basis functions in all integration points
      double* ed;
      int index_e = ref_map_shapeset.get_edge_index(k,0,ii+2);
      ref_map_pss.set_active_shape(index_e);
      ref_map_pss.set_quad_order(mo2);
      ed = ref_map_pss.get_fn_values();

      for (int m = 0; m < 2; m++)  //part 0 or 1
        for (int j = 0; j < np; j++)
          old[m][j] += proj[e->nvert + k * (order-1) + ii][m] * ed[j];
    }
  }
}


static void calc_bubble_projection(Element* e, Nurbs** nurbs, int order, double2* proj)
{
  ref_map_pss.set_active_element(e);

  int i, j, k;
  int mo2 = quad2d.get_max_order();
  int np = quad2d.get_num_points(mo2);
  int qo = e->is_quad() ? H2D_MAKE_QUAD_ORDER(order, order) : order;
  int nb = ref_map_shapeset.get_num_bubbles(qo);

  AUTOLA_OR(double2, fn, np);
  memset(fn, 0, sizeof(double2) * np);

  double* rhside[2];
  double* old[2];
  for (i = 0; i < 2; i++) {
    rhside[i] = new double[nb];
    old[i] = new double[np];
    memset(rhside[i], 0, sizeof(double) * nb);
    memset(old[i], 0, sizeof(double) * np);
  }

  // compute known part of projection (vertex and edge part)
  old_projection(e, order, proj, old);

  // fn values of both components of nonpolynomial function
  double3* pt = quad2d.get_points(mo2);
  for (j = 0; j < np; j++)  // over all integration points
  {
    double2 a;
    a[0] = ctm.m[0] * pt[j][0] + ctm.t[0];
    a[1] = ctm.m[1] * pt[j][1] + ctm.t[1];
    calc_ref_map(e, nurbs, a[0], a[1], fn[j]);
  }

  double2* result = proj + e->nvert + e->nvert * (order - 1);
  for (k = 0; k < 2; k++)
  {
    for (i = 0; i < nb; i++) // loop over bubble basis functions
    {
      // bubble basis functions in all integration points
      double *bfn;
      int index_i = ref_map_shapeset.get_bubble_indices(qo)[i];
      ref_map_pss.set_active_shape(index_i);
      ref_map_pss.set_quad_order(mo2);
      bfn = ref_map_pss.get_fn_values();

      for (j = 0; j < np; j++) // over all integration points
        rhside[k][i] += pt[j][2] * (bfn[j] * (fn[j][k] - old[k][j]));
    }

    // solve
    if (e->nvert == 3)
      cholsl(bubble_proj_matrix_tri, nb, bubble_tri_p, rhside[k], rhside[k]);
    else
      cholsl(bubble_proj_matrix_quad, nb, bubble_quad_p, rhside[k], rhside[k]);

    for (i = 0; i < nb; i++)
      result[i][k] = rhside[k][i];
  }

  for (i = 0; i < 2; i++) {
    delete [] rhside[i];
    delete [] old[i];
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////

static void ref_map_projection(Element* e, Nurbs** nurbs, int order, double2* proj)
{
  // vertex part
  for (unsigned int i = 0; i < e->nvert; i++)
  {
    proj[i][0] = e->vn[i]->x;
    proj[i][1] = e->vn[i]->y;
  }

  if (e->cm->toplevel == false)
    e = e->cm->parent;

  // edge part
  for (int edge = 0; edge < (int)e->nvert; edge++)
    calc_edge_projection(e, edge, nurbs, order, proj);

  //bubble part
  calc_bubble_projection(e, nurbs, order, proj);
}


void CurvMap::update_refmap_coefs(Element* e)
{
  ref_map_pss.set_quad_2d(&quad2d);
  //ref_map_pss.set_active_element(e);

  // calculation of projection matrices
  if (edge_proj_matrix == NULL) precalculate_cholesky_projection_matrix_edge();
  if (bubble_proj_matrix_tri == NULL) precalculate_cholesky_projection_matrices_bubble();

  ref_map_pss.set_mode(e->get_mode());
  ref_map_shapeset.set_mode(e->get_mode());

  // allocate projection coefficients
  int nv = e->nvert;
  int ne = order - 1;
  int qo = e->is_quad() ? H2D_MAKE_QUAD_ORDER(order, order) : order;
  int nb = ref_map_shapeset.get_num_bubbles(qo);
  nc = nv + nv*ne + nb;
  if (coefs != NULL) delete [] coefs;
  coefs = new double2[nc];

  // WARNING: do not change the format of the array 'coefs'. If it changes,
  // RefMap::set_active_element() has to be changed too.

  Nurbs** nurbs;
  if (toplevel == false)
  {
    ref_map_pss.set_active_element(e);
    ref_map_pss.set_transform(part);
    nurbs = parent->cm->nurbs;
  }
  else
  {
    ref_map_pss.reset_transform();
    nurbs = e->cm->nurbs;
  }
  ctm = *(ref_map_pss.get_ctm());
  ref_map_pss.reset_transform(); // fixme - do we need this?

  // calculation of new projection coefficients
  ref_map_projection(e, nurbs, order, coefs);
}


void CurvMap::get_mid_edge_points(Element* e, double2* pt, int n)
{
  Nurbs** nurbs = this->nurbs;
  Transformable tran;
  tran.set_active_element(e);

  if (toplevel == false)
  {
    tran.set_transform(part);
    e = e->cm->parent;
    nurbs = e->cm->nurbs;
  }

  ctm = *(tran.get_ctm());
  double xi_1, xi_2;
  for (int i = 0; i < n; i++)
  {
    xi_1 = ctm.m[0] * pt[i][0] + ctm.t[0];
    xi_2 = ctm.m[1] * pt[i][1] + ctm.t[1];
    calc_ref_map(e, nurbs, xi_1, xi_2, pt[i]);
  }
}


void Nurbs::unref()
{
  if (!--ref) // fixme: possible leak, we need ~Nurbs too
  {
    delete [] pt;
    delete [] kv;
    delete this;
  }
}


CurvMap::CurvMap(CurvMap* cm)
{
  memcpy(this, cm, sizeof(CurvMap));
  coefs = new double2[nc];
  memcpy(coefs, cm->coefs, sizeof(double2) * nc);

  if (toplevel)
    for (int i = 0; i < 4; i++)
      if (nurbs[i] != NULL)
        nurbs[i]->ref++;
}


CurvMap::~CurvMap()
{
  if (coefs != NULL)
    delete [] coefs;

  if (toplevel)
    for (int i = 0; i < 4; i++)
      if (nurbs[i] != NULL)
        nurbs[i]->unref();
}
