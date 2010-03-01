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
#include "solution.h"
#include "linsystem.h"
#include "refmap.h"
#include "shapeset_h1_all.h"
#include "quad_all.h"
#include "integrals_h1.h"
#include "matrix.h"
#include "adapt_ortho_h1.h"
#include "traverse.h"
#include "norm.h"

using namespace std;

H1OrthoHP::H1OrthoHP(int num, ...)
{
  this->num = num;

  va_list ap;
  va_start(ap, num);
  for (int i = 0; i < num; i++)
    spaces[i] = va_arg(ap, Space*);
  va_end(ap);

  for (int i = 0; i < num; i++)
    for (int j = 0; j < num; j++)
    {
      if (i == j) {
        form[i][j] = h1_form<double, scalar>;
        ord[i][j]  = h1_form<Ord, Ord>;
      }
      else {
        form[i][j] = NULL;
        ord[i][j]  = NULL;
      }
    }

  memset(errors, 0, sizeof(errors));
  esort = NULL;
  have_errors = false;
}


H1OrthoHP::~H1OrthoHP()
{
  for (int i = 0; i < num; i++)
    if (errors[i] != NULL)
      delete [] errors[i];

  if (esort != NULL) delete [] esort;
}


//// orthonormal base construction /////////////////////////////////////////////////////////////////

double3** H1OrthoHP::obase[2][9];
int  H1OrthoHP::basecnt[2][11];
bool H1OrthoHP::obase_ready = false;

int H1OrthoHP::build_shape_inxs(const int mode, H1Shapeset& shapeset, int idx[121]) {
  shapeset.set_mode(mode);

  // obtain a list of all shape functions up to the order 10, from lowest to highest order
  int n = 0;
  int nv = mode ? 4 : 3;
  int num_sons = mode ? 8 : 4;
  for (int i = 0; i < nv; i++)
    idx[n++] = shapeset.get_vertex_index(i);
  basecnt[mode][0] = 0;
  basecnt[mode][1] = n;

  for (int i = 2; i <= 10; i++)
  {
    for (int j = 0; j < nv; j++)
      idx[n++] = shapeset.get_edge_index(j, 0, i);

    int ii = mode ? make_quad_order(i, i) : i;
    int nb = shapeset.get_num_bubbles(ii);
    int* bub = shapeset.get_bubble_indices(ii);
    for (int j = 0; j < nb; j++)
    {
      int o = shapeset.get_order(bub[j]);
      if (get_h_order(o) == i || get_v_order(o) == i)
        idx[n++] = bub[j];
    }
    basecnt[mode][i] = n;
  }

  return n;
}

void H1OrthoHP::calc_ortho_base()
{
  int i, j, k, l, m, np, r;
  int n, idx[121];

  H1Shapeset shapeset;

  // allocate the orthonormal base tables - these are simply the values of the
  // orthonormal functions in integration points; we store the basic functions
  // plus four son cut-outs of them (i.e. 5 times)
  for (i = 0; i < 9; i++)
  {
    if ((i < 4) || (i >= 8))
      obase[0][i] = new_matrix<double3>(66, 79); // tri
    obase[1][i] = new_matrix<double3>(121, 121); // quad
  }

  // repeat for triangles and quads
  for (m = 0; m <= 1; m++)
  {
    //build indices
    n = build_shape_inxs(m, shapeset, idx);

    // obtain their values for integration rule 20
    g_quad_2d_std.set_mode(m);
    np = g_quad_2d_std.get_num_points(20);
    double3* pt = g_quad_2d_std.get_points(20);

    for (i = 0; i < n; i++)
      for (j = 0; j < np; j++)
        for (k = 0; k < 3; k++)
          obase[m][8][i][j][k] = shapeset.get_value(k, idx[i], pt[j][0], pt[j][1], 0);

    int num_sons = m ? 8 : 4;
    for (l = 0; l < num_sons; l++)
    {
      Trf* tr = (m ? quad_trf : tri_trf) + l;
      for (i = 0; i < n; i++)
        for (j = 0; j < np; j++)
        {
          double x = tr->m[0]*pt[j][0] + tr->t[0],
                 y = tr->m[1]*pt[j][1] + tr->t[1];
          for (k = 0; k < 3; k++)
            obase[m][l][i][j][k] = shapeset.get_value(k, idx[i], x, y, 0);
        }
    }

    // orthonormalize the basis functions
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < i; j++)
      {
        double prod = 0.0;
        for (k = 0; k < np; k++) {
          double sum = 0.0;
          for (r = 0; r < 3; r++)
            sum += obase[m][8][i][k][r] * obase[m][8][j][k][r];
          prod += pt[k][2] * sum;
        }

        for (l = 0; l < 9; l++)
          if (m || l < 4 || l >= 8)
            for (k = 0; k < np; k++)
              for (r = 0; r < 3; r++)
                obase[m][l][i][k][r] -= prod * obase[m][l][j][k][r];
      }

      double norm = 0.0;
      for (k = 0; k < np; k++) {
        double sum = 0.0;
        for (r = 0; r < 3; r++)
          sum += sqr(obase[m][8][i][k][r]);
        norm += pt[k][2] * sum;
      }
      norm = sqrt(norm);

      for (l = 0; l < 9; l++)
        if (m || l < 4 || l >= 8)
          for (k = 0; k < np; k++)
            for (r = 0; r < 3; r++)
              obase[m][l][i][k][r] /= norm;
    }

    // check the orthonormal base
/*    if (m) {
    for (i = 0; i < n; i++)
      for (j = 0; j < n; j++)
      {
        double check = 0.0;
        for(int son = 4; son < 6; son++ )
          for (k = 0; k < np; k++)
            check += pt[k][2] * (obase[m][son][i][k][0] * obase[m][son][j][k][0] +
                                 obase[m][son][i][k][1] * obase[m][son][j][k][1] +
                                 obase[m][son][i][k][2] * obase[m][son][j][k][2]);
        check *= 0.5;
        if ((i == j && fabs(check - 1.0) > 1e-8) || (i != j && fabs(check) > 1e-8))
          warn("Not orthonormal: base %d times base %d = %g", i, j , check);
      }
    }*/
  }
  obase_ready = true;
}


void H1OrthoHP::free_ortho_base()
{
  if (!obase_ready) return;

  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 2; j++)
      delete [] obase[j][i];

  obase_ready = false;
}


//// optimal refinement search /////////////////////////////////////////////////////////////////////

void H1OrthoHP::calc_projection_errors(Element* e, int order, Solution* rsln,
                                       double herr[8][11], double perr[11])
{
  int i, j, s, k, r, son;
  int m = e->get_mode();
  double error;
  scalar prod;

  if (!obase_ready) calc_ortho_base();

  // select quadrature, obtain integration points and weights
  Quad2D* quad = &g_quad_2d_std;
  quad->set_mode(m);
  rsln->set_quad_2d(quad);
  double3* pt = quad->get_points(20);
  int np = quad->get_num_points(20);

  // everything is done on the reference domain
  // -- no reference mapping, no transformations
  rsln->enable_transform(false);

  // obtain reference solution values on all four refined sons
  scalar* rval[4][3];
  Element* base = rsln->get_mesh()->get_element(e->id);
  assert(!base->active);
  for (son = 0; son < 4; son++)
  {
    Element* e = base->sons[son];
    assert(e != NULL);
    rsln->set_active_element(e);
    rsln->set_quad_order(20);
    rval[son][0] = rsln->get_fn_values();
    rval[son][1] = rsln->get_dx_values();
    rval[son][2] = rsln->get_dy_values();
  }

  // h-cadidates: calculate products of the reference solution with orthonormal basis
  // functions on son elements, obtaining (partial) projections and their errors
  // the error is scaled by 4 (error of four sons of an element of summed togethers; error for every element is evaluted in a reference domain
  scalar3 proj[4][121];
  for (son = 0; son < 4; son++)
  {
    memset(proj[0], 0, sizeof(proj[0]));
    for (i = 1; i <= order; i++)
    {
      // update the projection to the current order
      for (j = basecnt[m][i-1]; j < basecnt[m][i]; j++)
      {
        for (k = 0, prod = 0.0; k < np; k++)
          prod += pt[k][2] * (rval[son][0][k] * obase[m][8][j][k][0] +
                              rval[son][1][k] * obase[m][8][j][k][1] +
                              rval[son][2][k] * obase[m][8][j][k][2]);

        for (k = 0; k < np; k++)
          for (r = 0; r < 3; r++)
            proj[0][k][r] += obase[m][8][j][k][r] * prod;
      }

      // calculate the error of the projection
      for (k = 0, error = 0.0; k < np; k++)
        error += pt[k][2] * (sqr(rval[son][0][k] - proj[0][k][0]) +
                             sqr(rval[son][1][k] - proj[0][k][1]) +
                             sqr(rval[son][2][k] - proj[0][k][2]));
      herr[son][i] = error;
    }
  }

  // aniso-candidates: calculate projections and their errors (only quadrilaterals)
  // the error is not scaled
  if (m)
  {
    const double mx[4] = { 2.0, 2.0, 1.0, 1.0};
    const double my[4] = { 1.0, 1.0, 2.0, 2.0};
    const int sons[4][2] = { {0,1}, {3,2}, {0,3}, {1,2} };
    const int tr[4][2]   = { {6,7}, {6,7}, {4,5}, {4,5} };

    for (son = 0; son < 4; son++) // 2 sons for vertical split, 2 sons for horizontal split
    {
      memset(proj, 0, sizeof(proj));
      for (i = 1; i <= order+1; i++)  // h-candidates: max order equals to original element order+1
      {
        // update the projection to the current order
        for (j = basecnt[m][i-1]; j < basecnt[m][i]; j++)
        {
          scalar prod = 0.0;
          for (s = 0; s < 2; s++) // each son has 2 subsons (regular square sons)
            for (k = 0; k < np; k++)
              prod += pt[k][2] * ( rval[sons[son][s]][0][k]           * obase[m][tr[son][s]][j][k][0] +
                                   rval[sons[son][s]][1][k] * mx[son] * obase[m][tr[son][s]][j][k][1] +
                                   rval[sons[son][s]][2][k] * my[son] * obase[m][tr[son][s]][j][k][2]);
          prod *= 0.5; //compensate the fact that values are a sum of itegral of two elements in a ref. domain

          for (s = 0; s < 2; s++)
            for (k = 0; k < np; k++)
              for (r = 0; r < 3; r++)
                proj[s][k][r] += prod * obase[m][tr[son][s]][j][k][r];
        }

        // calculate the error of the projection
        for (s = 0, error = 0.0; s < 2; s++)
          for (k = 0; k < np; k++)
            error += pt[k][2] * (sqr(rval[sons[son][s]][0][k]           - proj[s][k][0]) +
                                 sqr(rval[sons[son][s]][1][k] * mx[son] - proj[s][k][1]) +
                                 sqr(rval[sons[son][s]][2][k] * my[son] - proj[s][k][2]));
        herr[4 + son][i] = error * 0.5; //compensate the fact that values are a sum of itegral of two elements in a ref. domain
      }
    }
  }

  // p-candidates: calculate projections and their errors
  memset(proj, 0, sizeof(proj));
  for (i = 1; i <= std::min(order+2, 10); i++)
  {
    // update the projection to the current order
    for (j = basecnt[m][i-1]; j < basecnt[m][i]; j++)
    {
      scalar prod = 0.0;
      for (son = 0; son < 4; son++)
      {
        // (transforming to the quarter of the reference element)
        double mm = (e->is_triangle() && son == 3) ? -2.0 : 2.0;

        for (k = 0; k < np; k++)
        {
          prod += pt[k][2] * (rval[son][0][k] *      obase[m][son][j][k][0] +
                              rval[son][1][k] * mm * obase[m][son][j][k][1] +
                              rval[son][2][k] * mm * obase[m][son][j][k][2]);
        }
      }
      prod *= 0.25; //compensate the fact that values are a sum of itegral of four elements in a ref. domain

      for (son = 0; son < 4; son++)
        for (k = 0; k < np; k++)
          for (r = 0; r < 3; r++)
            proj[son][k][r] += prod * obase[m][son][j][k][r];
    }

    // calculate the error of the projection
    for (son = 0, error = 0.0; son < 4; son++)
    {
      double mm = (e->is_triangle() && son == 3) ? -2.0 : 2.0;

      for (k = 0; k < np; k++)
        error += pt[k][2] * (sqr(rval[son][0][k]      - proj[son][k][0]) +
                             sqr(rval[son][1][k] * mm - proj[son][k][1]) +
                             sqr(rval[son][2][k] * mm - proj[son][k][2]));
    }
    perr[i] = error * 0.25; //compensate the fact that values are a sum of itegral of four elements in a ref. domain
  }

  rsln->enable_transform(true);

}

H1OrthoHP::Cand* H1OrthoHP::create_candidates(Element* e, int order, bool h_only, bool iso_only, int max_order, int* num_cand) {
  int n = 0;
  const int maxcand = 300;

  bool tri = e->is_triangle();

  // calculate default maximal order of elements
  // linear elements = 9
  // curvilinear elements = depends on iro_cache (how curved they are)
  if (max_order == -1)
    max_order = (20 - e->iro_cache)/2 - 1; // default
  else
    max_order = std::min( max_order, (20 - e->iro_cache)/2 - 1); // user specified

  int min_order = 1;

  Cand* cand = new Cand[maxcand];

  #define make_p_cand(q) { \
    assert(n < maxcand);   \
    cand[n].split = -1; \
    cand[n].p[1] = cand[n].p[2] = cand[n].p[3] = 0; \
    cand[n++].p[0] = (q); }

  #define make_hp_cand(q0, q1, q2, q3) { \
    assert(n < maxcand);  \
    cand[n].split = 0; \
    cand[n].p[0] = (q0); \
    cand[n].p[1] = (q1); \
    cand[n].p[2] = (q2); \
    cand[n++].p[3] = (q3); }

  #define make_ani_cand(q0, q1, iso) { \
    assert(n < maxcand);  \
    cand[n].split = iso; \
    cand[n].p[2] = cand[n].p[3] = 0; \
    cand[n].p[0] = (q0); \
    cand[n++].p[1] = (q1); }

  if (h_only)
  {
    make_p_cand(order);
    make_hp_cand(order, order, order, order);
    if ((!tri) && (e->iro_cache < 8) && !iso_only) {
      make_ani_cand(order, order, 1);
      make_ani_cand(order, order, 2);
    }
  }
  else {
    // prepare p-candidates
    int p0, p1 = std::min(max_order, order+2);
    for (p0 = order; p0 <= p1; p0++) {
      make_p_cand(p0);
    }

    //prepare hp-candidates
    p0 = std::max(min_order, (order+1) / 2);
    p1 = std::max(min_order, std::min(p0 + 3, order));
    int q0, q1, q2, q3;
    for (q0 = p0; q0 <= p1; q0++)
      for (q1 = p0; q1 <= p1; q1++)
        for (q2 = p0; q2 <= p1; q2++)
          for (q3 = p0; q3 <= p1; q3++)
            make_hp_cand(q0, q1, q2, q3);

    //prepare anisotropic candidates
    //only for quadrilaterals
    //too distorted (curved) elements cannot have aniso refinement (produces even worse elements)
    if ((!tri) && (e->iro_cache < 8) && !iso_only) {
      p0 = 2 * (order+1) / 3;
      int p_max = std::min(max_order, order+1);
      p1 = std::min(p0 + 3, p_max);
      for (q0 = p0; q0 <= p1; q0++)
        for (q1 = p0; q1 <= p1; q1++) {
          if ((q0 < order+1) || (q1 < order+1)) {
            make_ani_cand(q0, q1, 1);
            make_ani_cand(q0, q1, 2);
          }
        }
    }
  }

  *num_cand = n;
  return cand;
}

int H1OrthoHP::evalute_candidates(Cand* cand, int num_cand, Element* e, int order, Solution* rsln, double* avg_error, double* dev_error) {
  bool tri = e->is_triangle();

  // calculate (partial) projection errors
  double herr[8][11], perr[11];
  calc_projection_errors(e, order, rsln, herr, perr);

  //evaluate errors and dofs
  double sum_err = 0.0;
  double sum_sqr_err = 0.0;
  int num_processed = 0;
  for (int i = 0; i < num_cand; i++)
  {
    Cand* c = cand + i;
    if (c->split == 0)
    {
      c->error = 0.0;
      c->dofs = tri ? 6 : 9;
      for (int j = 0; j < 4; j++)
      {
        int o = c->p[j];
        c->error += herr[j][o]; // * 0.25; // spravny vypocet chyby, ??? candidate error is composed of four sons
        if (tri) {
          c->dofs += (o-2)*(o-1)/2;
          if (j < 3) c->dofs += std::min(o, c->p[3])-1 + 2*(o-1);
        }
        else {
          c->dofs += sqr(o)-1;
          c->dofs += /* 2 * */ std::min(o, c->p[j>0 ? j-1 : 3]) - 1; //??? should it be commented out or not?
        }
      }
    }
    else if (c->split == 1 || c->split == 2)  // aniso splits
    {
      c->dofs  = 6 /* vertex */ + 3*(c->p[0] - 1 + c->p[1] - 1); // edge fns
      c->dofs += std::min(c->p[0], c->p[1]) - 1; // common edge
      c->dofs += sqr(c->p[0] - 1) + sqr(c->p[1] - 1); // bubbles
      c->error = 0.0;
      for (int j = 0; j < 2; j++)
        c->error += herr[(c->split == 1) ? j+4 : j+6][c->p[j]]; // * 0.5;  // spravny vypocet chyby, ??? average of errors on splot element (sons)

    }
    else
    {
      int o = c->p[0];
      c->error = perr[o];
      c->dofs  = tri ? (o+1)*(o+2)/2 : sqr(o+1);
    }
    c->error = sqrt(c->error);

    if (!i || c->error <= cand[0].error)
    {
      sum_err += log(c->error);
      sum_sqr_err += sqr(log(c->error));
      num_processed++;
    }
  }
 
  *avg_error = sum_err / num_processed;  // mean
  *dev_error = sqrt(sum_sqr_err/num_processed - sqr(*avg_error)); // deviation is square root of variance
  return num_processed;
}

void H1OrthoHP::select_best_candidate(const Cand* cand, const int num_cand, Element* e, const double avg_error, const double dev_error, int* selected_cand, int* selected_h_cand, double conv_exp) {
  // select an above-average candidate with the steepest error decrease
  int imax = 0, h_imax = 0;
  double score, maxscore = 0.0, h_maxscore = 0.0;
  for (int i = 1; i < num_cand; i++)
  {
    if ((log(cand[i].error) < (avg_error + dev_error)) && (cand[i].dofs > cand[0].dofs))
    {
      score = (log(cand[0].error) - log(cand[i].error)) / 
	       //(pow(cand[i].dofs, conv_exp) - pow(cand[0].dofs, conv_exp));
               pow(cand[i].dofs - cand[0].dofs, conv_exp);

      if (score > maxscore) { maxscore = score; imax = i; }
      if ((cand[i].split == 0) && (score > h_maxscore)) { h_maxscore = score; h_imax = i; }
    }
  }

  *selected_cand = imax;
  *selected_h_cand = h_imax;
}

int H1OrthoHP::get_optimal_refinement(Element* e, int order, Solution* rsln, int& split, int4 p, int4 q,
				      bool h_only, bool iso_only, double conv_exp, int max_order)
{
  //decode order
  order = std::max(get_h_order(order), get_v_order(order));

  //build candidates
  int num_cand = -1;
  Cand* cand = create_candidates(e, order, h_only, iso_only, max_order, &num_cand);

  // evaluate candidates (sum partial projection errors, calculate dofs)
  double avg_error, dev_error;
  evalute_candidates(cand, num_cand, e, order, rsln, &avg_error, &dev_error);

  //select candidate
  int inx_cand, inx_h_cand;
  select_best_candidate(cand, num_cand, e, avg_error, dev_error, &inx_cand, &inx_h_cand, conv_exp);

  //copy result to output
  split = cand[inx_cand].split;
  memcpy(p, cand[inx_cand].p, 4*sizeof(int));
  memcpy(q, cand[inx_h_cand].p, 4*sizeof(int));

  //clenaup
  delete[] cand;

  return inx_cand;
}


//// adapt /////////////////////////////////////////////////////////////////////////////////////////

bool H1OrthoHP::adapt(double thr, int strat, int adapt_type, bool iso_only, int regularize,
                      double conv_exp, int max_order, bool same_orders, double to_be_processed)
{
  if (!have_errors)
    error("Element errors have to be calculated first, see calc_error().");

  int i, j, l;
  int max_id = -1;
  Mesh* meshes[10];
  for (j = 0; j < num; j++) {
    meshes[j] = spaces[j]->get_mesh();
    rsln[j]->set_quad_2d(&g_quad_2d_std);
    rsln[j]->enable_transform(false);
    if (meshes[j]->get_max_element_id() > max_id)
      max_id = meshes[j]->get_max_element_id();
  }

  AUTOLA2_OR(int, idx, max_id + 1, num + 1);
  for(j = 0; j < max_id; j++)
    for(l = 0; l < num; l++)
      idx[j][l] = -1; // element not refined

  //int nref = nact;
  double err0 = 1000.0;
  double processed_error = 0.0;
  bool h_only = adapt_type == 1 ? true : false;

  vector<ElementToRefine> elem_inx_to_proc; //list of indices of elements that are going to be processed
  elem_inx_to_proc.reserve(nact);

  //adaptivity loop
  double error_threshod = -1; //an error threshold that breaks the adaptivity loop in a case of strategy 1
  int num_exam_elem = 0; //a number of examined elements
  int num_ignored_elem = 0; //a number of ignored elements
  int num_not_changed = 0; //a number of element that were not changed
  int num_priority_elem = 0; //a number of elements that were processed using priority queue

  int inx_regular_element = 0;
  while (inx_regular_element < nact || !priority_esort.empty())
  {
    int id, comp, inx_element;

    //get element identification
    if (priority_esort.empty()) {
      id = esort[inx_regular_element].id;
      comp = esort[inx_regular_element].comp;
      inx_element = inx_regular_element;
      inx_regular_element++;
    }
    else {
      id = priority_esort.front().id;
      comp = priority_esort.front().comp;
      inx_element = -1;
      priority_esort.pop();
      num_priority_elem++;
    }
    num_exam_elem++;

    //get info linked with the element
    double err = errors[comp][id];
    Mesh* mesh = meshes[comp];
    Element* e = mesh->get_element(id);

    if (!ignore_element_adapt(inx_element, mesh, e)) {

      //use error of the first refined element to calculate the threshold for strategy 1
      if (elem_inx_to_proc.empty())
        error_threshod = thr * err;

      //stop the loop only if processing regular elements
      if (inx_element >= 0) {

        // first refinement strategy:
        // refine elements until prescribed amount of error is processed
        // if more elements have similar error refine all to keep the mesh symmetric
        if ((strat == 0) && (processed_error > sqrt(thr) * total_err) && fabs((err - err0)/err0) > 1e-3) break;

        // second refinement strategy:
        // refine all elements whose error is bigger than some portion of maximal error
        if ((strat == 1) && (err < error_threshod)) break;

        if ((strat == 2) && (err < thr)) break;

        if ((strat == 3) &&
          ( (err < error_threshod) ||
          ( processed_error > 1.5 * to_be_processed )) ) break;
      }

      // p-adaptivity
      ElementToRefine elem_ref(id, comp);
      int current = spaces[comp]->get_element_order(id);
      bool refined = false;
      if (adapt_type == 2) {
        elem_ref.split = -1;
        elem_ref.p[0] = elem_ref.q[0] = std::min(9, get_h_order(current) + 1);
        if (get_h_order(current) < elem_ref.p[0]) refined = true;
      }
      // h-adaptivity
      else if (adapt_type == 1 && iso_only) {
        elem_ref.split = 0;
        elem_ref.p[0] = elem_ref.p[1] = elem_ref.p[2] = elem_ref.p[3] = current;
        elem_ref.q[0] = elem_ref.q[1] = elem_ref.q[2] = elem_ref.q[3] = current;
        refined = true;
      }
      // hp-adaptivity
      else {
        int inx_candidate = get_optimal_refinement(e, current, rsln[comp], 
                                                   elem_ref.split, elem_ref.p, elem_ref.q, 
                                                   h_only, iso_only, conv_exp, max_order);
        if (inx_candidate != 0)
          refined = true;
      }

      //add to a list of elements that are going to be refined
      if (refined && can_adapt_element(mesh, e, elem_ref.split, elem_ref.p, elem_ref.q) ) {
        idx[id][comp] = (int)elem_inx_to_proc.size();
        elem_inx_to_proc.push_back(elem_ref);
        err0 = err;
        processed_error += err;
      }
      else
        num_not_changed++;
    }
    else {
      num_ignored_elem++;
    }
  }

  debug_log("I examined elements: %d", num_exam_elem);
  debug_log("  elements taken from priority queue: %d", num_priority_elem);
  debug_log("  ignored elements: %d", num_ignored_elem);
  debug_log("  not changed elements: %d", num_not_changed);
  debug_log("  elements to process: %d", elem_inx_to_proc.size());
  bool done = false;
  if (num_exam_elem == 0)
    done = true;
  else if (elem_inx_to_proc.empty())
  {
    warn("W none of the elements selected for refinement could be refined. Adaptivity step not successful, returning 'true'.\n");
    done = true;
  }

  int num_elem_to_proc = elem_inx_to_proc.size();
  for(int inx = 0; inx < num_elem_to_proc; inx++) {
    ElementToRefine& elem_ref = elem_inx_to_proc[inx];
    int current = get_h_order(spaces[elem_ref.comp]->get_element_order(elem_ref.id));
    
    int max_ref = elem_ref.split; // how all the other elements will be refined
    for (j = 0; j < num; j++)
    {
      if (max_ref == 0) break; // iso refinement is max what can be recieved
      if ((j != elem_ref.comp) && (meshes[j] == meshes[elem_ref.comp])) // components share the mesh
      {
        int ii = idx[elem_ref.id][j];
        if (ii >= 0) {
          const ElementToRefine& elem_ref_ii = elem_inx_to_proc[ii];
          if ((elem_ref_ii.split != max_ref) && (elem_ref_ii.split >= 0)) { // ii element refined, refinement differs from max_ref, ii element split
            if (((elem_ref_ii.split == 1) || (elem_ref_ii.split == 2)) && (max_ref == -1)) // the only case when aniso refinement
              max_ref = elem_ref_ii.split;
            else // otherwise isotropic refinement
              max_ref = 0;
          }
        }
      }
    }

    if (max_ref >= 0)
    {
      for (j = 0; j < num; j++)
      {
        if ((j != elem_ref.comp) && (meshes[j] == meshes[elem_ref.comp])) // components share the mesh
        {
          // change appropriately original element
          if (elem_ref.split != max_ref)
          {
            elem_ref.split = max_ref;
            if (elem_ref.split == 0)
              memcpy(elem_ref.p, elem_ref.q, 4*sizeof(int));
            else { // aniso refinements
              elem_ref.p[0] = h_only ? current : std::max(1, 2*(current+1)/3);
              elem_ref.p[1] = h_only ? current : std::max(1, 2*(current+1)/3);
            }
          }
          int ii = idx[elem_ref.id][j];
          current = get_h_order(spaces[j]->get_element_order(elem_ref.id));
          if (ii >= 0)
          {
            ElementToRefine& elem_ref_ii = elem_inx_to_proc[ii];
            if (elem_ref_ii.split != max_ref)
            {
              elem_ref_ii.split = max_ref;
              if (elem_ref_ii.split == 0)
                memcpy(elem_ref_ii.p, elem_ref_ii.q, 4*sizeof(int));
              else { // aniso refinements
                elem_ref_ii.p[0] = h_only ? current : std::max(1, 2*(current+1)/3);
                elem_ref_ii.p[1] = h_only ? current : std::max(1, 2*(current+1)/3);
              }
            }
          }
          if (ii < 0) // element not refined at all
          {
            ElementToRefine elem_ref_new(elem_ref.id, j);
            elem_ref_new.split = max_ref;
            if (elem_ref_new.split == 0)
              for (int r = 0; r < 4; r++)
                elem_ref_new.p[r] = h_only ? current : std::max(1, (current+1)/2);
            else { // aniso refinements
              elem_ref_new.p[0] = h_only ? current : std::max(1, 2*(current+1)/3);
              elem_ref_new.p[1] = h_only ? current : std::max(1, 2*(current+1)/3);
            }
            elem_inx_to_proc.push_back(elem_ref_new);
          }
        }
      }
    }
  }

  //apply refinements
  apply_refinements(meshes, &elem_inx_to_proc);

  if (same_orders)
  {
    Element* e;
    for (i = 0; i < num; i++)
    {
      for_all_active_elements(e, meshes[i])
      {
        int current = get_h_order(spaces[i]->get_element_order(e->id));
        for (j = 0; j < num; j++)
          if ((j != i) && (meshes[j] == meshes[i])) // components share the mesh
          {
            int o = get_h_order(spaces[j]->get_element_order(e->id));
            if (o > current) current = o;
          }
        spaces[i]->set_element_order(e->id, current);
      }
    }
  }

  // mesh regularization
  if (regularize >= 0)
  {
    if (regularize == 0)
    {
      regularize = 1;
      warn("Total mesh regularization is not supported in adaptivity. 1-irregular mesh is used instead.");
    }
    for (i = 0; i < num; i++)
    {
      int* parents;
      parents = meshes[i]->regularize(regularize);
      spaces[i]->distribute_orders(meshes[i], parents);
      delete [] parents;
    }
  }

  for (j = 0; j < num; j++)
    rsln[j]->enable_transform(true);


  verbose("I refined elements: %d", elem_inx_to_proc.size());
  have_errors = false;
  if (strat == 2 && done == true) have_errors = true; // space without changes

  return done;
}

void H1OrthoHP::apply_refinements(Mesh** meshes, std::vector<ElementToRefine>* elems_to_refine)
{
  for (vector<ElementToRefine>::const_iterator elem_ref = elems_to_refine->begin(); elem_ref != elems_to_refine->end(); elem_ref++) // go over elements to be refined
  {
    Element* e;
    e = meshes[elem_ref->comp]->get_element(elem_ref->id);

    if (elem_ref->split < 0)
      spaces[elem_ref->comp]->set_element_order(elem_ref->id, elem_ref->p[0]);
    else if (elem_ref->split == 0) {
      if (e->active)
        meshes[elem_ref->comp]->refine_element(elem_ref->id);
      for (int j = 0; j < 4; j++)
        spaces[elem_ref->comp]->set_element_order(e->sons[j]->id, elem_ref->p[j]);
    }
    else {
      if (e->active)
        meshes[elem_ref->comp]->refine_element(elem_ref->id, elem_ref->split);
      for (int j = 0; j < 2; j++)
        spaces[elem_ref->comp]->set_element_order(e->sons[ (elem_ref->split == 1) ? j : j+2 ]->id, elem_ref->p[j]);
    }
  }
}


///// Unrefinements /////////////////////////////////////////////////////////////////////////////////

void H1OrthoHP::unrefine(double thr)
{

  if (!have_errors)
    error("Element errors have to be calculated first, see calc_error().");

  Mesh* mesh[2];
  mesh[0] = spaces[0]->get_mesh();
  mesh[1] = spaces[1]->get_mesh();


  int k = 0;
  if (mesh[0] == mesh[1]) // single mesh
  {
    Element* e;
    for_all_inactive_elements(e, mesh[0])
    {
      bool found = true;
      for (int i = 0; i < 4; i++)
        if (e->sons[i] != NULL && ((!e->sons[i]->active) || (e->sons[i]->is_curved())))
      { found = false;  break; }

      if (found)
      {
        double sum1 = 0.0, sum2 = 0.0;
        int max1 = 0, max2 = 0;
        for (int i = 0; i < 4; i++)
          if (e->sons[i] != NULL)
        {
          sum1 += errors[0][e->sons[i]->id];
          sum2 += errors[1][e->sons[i]->id];
          int oo = spaces[0]->get_element_order(e->sons[i]->id);
          if (oo > max1) max1 = oo;
          oo = spaces[1]->get_element_order(e->sons[i]->id);
          if (oo > max2) max2 = oo;
        }
        if ((sum1 < thr * errors[esort[0].comp][esort[0].id]) &&
             (sum2 < thr * errors[esort[0].comp][esort[0].id]))
        {
          mesh[0]->unrefine_element(e->id);
          mesh[1]->unrefine_element(e->id);
          errors[0][e->id] = sum1;
          errors[1][e->id] = sum2;
          spaces[0]->set_element_order(e->id, max1);
          spaces[1]->set_element_order(e->id, max2);
          k++; // number of unrefined elements
        }
      }
    }
    for_all_active_elements(e, mesh[0])
    {
      for (int i = 0; i < 2; i++)
        if (errors[i][e->id] < thr/4 * errors[esort[0].comp][esort[0].id])
      {
        int oo = get_h_order(spaces[i]->get_element_order(e->id));
        spaces[i]->set_element_order(e->id, std::max(oo - 1, 1));
        k++;
      }
    }
  }
  else // multimesh
  {
    for (int m = 0; m < 2; m++)
    {
      Element* e;
      for_all_inactive_elements(e, mesh[m])
      {
        bool found = true;
        for (int i = 0; i < 4; i++)
          if (e->sons[i] != NULL && ((!e->sons[i]->active) || (e->sons[i]->is_curved())))
        { found = false;  break; }

        if (found)
        {
          double sum = 0.0;
          int max = 0;
          for (int i = 0; i < 4; i++)
            if (e->sons[i] != NULL)
          {
            sum += errors[m][e->sons[i]->id];
            int oo = spaces[m]->get_element_order(e->sons[i]->id);
            if (oo > max) max = oo;
          }
          if ((sum < thr * errors[esort[0].comp][esort[0].id]))
          //if ((sum < 0.1 * thr))
          {
            mesh[m]->unrefine_element(e->id);
            errors[m][e->id] = sum;
            spaces[m]->set_element_order(e->id, max);
            k++; // number of unrefined elements
          }
        }
      }
      for_all_active_elements(e, mesh[m])
      {
        if (errors[m][e->id] < thr/4 * errors[esort[0].comp][esort[0].id])
        {
          int oo = get_h_order(spaces[m]->get_element_order(e->id));
          spaces[m]->set_element_order(e->id, std::max(oo - 1, 1));
          k++;
        }
      }
    }
  }
  verbose("Unrefined %d elements.", k);
  have_errors = false;
}

//// error calculation /////////////////////////////////////////////////////////////////////////////

double** H1OrthoHP::cmp_err;
int H1OrthoHP::compare(const void* p1, const void* p2) {
  const ElementReference& e1 = *((const ElementReference*)p1);
  const ElementReference& e2 = *((const ElementReference*)p2);
  return cmp_err[e1.comp][e1.id] < cmp_err[e2.comp][e2.id] ? 1 : -1;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void H1OrthoHP::set_biform(int i, int j, biform_val_t bi_form, biform_ord_t bi_ord)
{
  if (i < 0 || i >= num || j < 0 || j >= num)
    error("Invalid equation number.");

  form[i][j] = bi_form;
  ord[i][j] = bi_ord;
}


scalar H1OrthoHP::eval_error(biform_val_t bi_fn, biform_ord_t bi_ord,
                             MeshFunction *sln1, MeshFunction *sln2, MeshFunction *rsln1, MeshFunction *rsln2,
                             RefMap *rv1,        RefMap *rv2,        RefMap *rrv1,        RefMap *rrv2)
{
  // determine the integration order
  int inc = (rsln1->get_num_components() == 2) ? 1 : 0;
  Func<Ord>* ou = init_fn_ord(rsln1->get_fn_order() + inc);
  Func<Ord>* ov = init_fn_ord(rsln2->get_fn_order() + inc);

  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  Ord o = bi_ord(1, &fake_wt, ou, ov, fake_e, NULL);
  int order = rrv1->get_inv_ref_order();
  order += o.get_order();
  limit_order(order);

  ou->free_ord(); delete ou;
  ov->free_ord(); delete ov;
  delete fake_e;

  // eval the form
  Quad2D* quad = sln1->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // init geometry and jacobian*weights
  Geom<double>* e = init_geom_vol(rrv1, order);
  double* jac = rrv1->get_jacobian(order);
  double* jwt = new double[np];
  for(int i = 0; i < np; i++)
    jwt[i] = pt[i][2] * jac[i];

  // function values and values of external functions
  Func<scalar>* err1 = init_fn(sln1, rv1, order);
  Func<scalar>* err2 = init_fn(sln2, rv2, order);
  Func<scalar>* v1 = init_fn(rsln1, rrv1, order);
  Func<scalar>* v2 = init_fn(rsln2, rrv2, order);

  for (int i = 0; i < np; i++)
  {
    err1->val[i] = err1->val[i] - v1->val[i];
    err1->dx[i] = err1->dx[i] - v1->dx[i];
    err1->dy[i] = err1->dy[i] - v1->dy[i];
    err2->val[i] = err2->val[i] - v2->val[i];
    err2->dx[i] = err2->dx[i] - v2->dx[i];
    err2->dy[i] = err2->dy[i] - v2->dy[i];
  }

  scalar res = bi_fn(np, jwt, err1, err2, e, NULL);

  e->free(); delete e;
  delete [] jwt;
  err1->free_fn(); delete err1;
  err2->free_fn(); delete err2;
  v1->free_fn(); delete v1;
  v2->free_fn(); delete v2;

  return res;
}


scalar H1OrthoHP::eval_norm(biform_val_t bi_fn, biform_ord_t bi_ord,
                            MeshFunction *rsln1, MeshFunction *rsln2, RefMap *rrv1, RefMap *rrv2)
{
  // determine the integration order
  int inc = (rsln1->get_num_components() == 2) ? 1 : 0;
  Func<Ord>* ou = init_fn_ord(rsln1->get_fn_order() + inc);
  Func<Ord>* ov = init_fn_ord(rsln2->get_fn_order() + inc);

  double fake_wt = 1.0;
  Geom<Ord>* fake_e = init_geom_ord();
  Ord o = bi_ord(1, &fake_wt, ou, ov, fake_e, NULL);
  int order = rrv1->get_inv_ref_order();
  order += o.get_order();
  limit_order(order);

  ou->free_ord(); delete ou;
  ov->free_ord(); delete ov;
  delete fake_e;

  // eval the form
  Quad2D* quad = rsln1->get_quad_2d();
  double3* pt = quad->get_points(order);
  int np = quad->get_num_points(order);

  // init geometry and jacobian*weights
  Geom<double>* e = init_geom_vol(rrv1, order);
  double* jac = rrv1->get_jacobian(order);
  double* jwt = new double[np];
  for(int i = 0; i < np; i++)
    jwt[i] = pt[i][2] * jac[i];

  // function values
  Func<scalar>* v1 = init_fn(rsln1, rrv1, order);
  Func<scalar>* v2 = init_fn(rsln2, rrv2, order);

  scalar res = bi_fn(np, jwt, v1, v2, e, NULL);

  e->free(); delete e;
  delete [] jwt;
  v1->free_fn(); delete v1;
  v2->free_fn(); delete v2;

  return res;
}


double H1OrthoHP::calc_error(MeshFunction* sln, MeshFunction* rsln)
{
  if (num != 1) error("Wrong number of solutions.");

  return calc_error_n(1, sln, rsln);
}


double H1OrthoHP::calc_error_2(MeshFunction* sln1, MeshFunction* sln2, MeshFunction* rsln1, MeshFunction* rsln2)
{
  if (num != 2) error("Wrong number of solutions.");

  return calc_error_n(2, sln1, sln2, rsln1, rsln2);
}


double H1OrthoHP::calc_error_n(int n, ...)
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
      RefMap* rmi = sln[i]->get_refmap();
      RefMap* rrmi = rsln[i]->get_refmap();
      for (j = 0; j < num; j++)
      {
        RefMap* rmj = sln[j]->get_refmap();
        RefMap* rrmj = rsln[j]->get_refmap();
        double e, t;
        if (form[i][j] != NULL)
        {
          #ifndef COMPLEX
          e = fabs(eval_error(form[i][j], ord[i][j], sln[i], sln[j], rsln[i], rsln[j], rmi, rmj, rrmi, rrmj));
          t = fabs(eval_norm(form[i][j], ord[i][j], rsln[i], rsln[j], rrmi, rrmj));
          #else
          e = std::abs(eval_error(form[i][j], ord[i][j], sln[i], sln[j], rsln[i], rsln[j], rmi, rmj, rrmi, rrmj));
          t = std::abs(eval_norm(form[i][j], ord[i][j], rsln[i], rsln[j], rrmi, rrmj));
          #endif

          norms[i] += t;
          total_norm  += t;
          total_error += e;
          errors[i][ee[i]->id] += e;
        }
      }
    }
  }
  trav.finish();

  //prepare an ordered list of elements according to an error
  sort_elements_by_error(meshes);

  have_errors = true;
  total_err = total_error/* / total_norm*/;
  return sqrt(total_error / total_norm);
}

void H1OrthoHP::sort_elements_by_error(Mesh** meshes) {
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
//       errors[i][e->id] /= norms[i];
// ??? needed or not ???
// when norms of 2 components are very different it can help (microwave heating)
// navier-stokes on different meshes work only without
    }

  //sort
  assert(inx == nact);
  cmp_err = errors;
  qsort(esort, nact, sizeof(ElementReference), compare);
}
