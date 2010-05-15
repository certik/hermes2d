#define H2D_REPORT_WARN
#define H2D_REPORT_INFO
#define H2D_REPORT_VERBOSE
#define H2D_REPORT_FILE "application.log"

#include "hermes2d.h"
#include "solver_umfpack.h"
#include "h1_adapt_norm_maxim.h"

double H1AdaptNormMaxim::norm_fn_inf(MeshFunction* sln, RefMap* ru) {
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

scalar H1AdaptNormMaxim::eval_norm(biform_val_t bi_fn, biform_ord_t bi_ord, MeshFunction *rsln1, MeshFunction *rsln2, RefMap *rrv1, RefMap *rrv2) {
  error_if(rsln1 != rsln2, "Multiple components not supported.");
  error_if(rrv1 != rrv2, "Multiple components not supported.");

  switch (norm_type) {
    case NORM_EUCLEDIAN:
      return H1Adapt::eval_norm(bi_fn, bi_ord, rsln1, rsln2, rrv1, rrv2);
    case NORM_MAXIMA:
      return norm_fn_inf(rsln1, rrv1);
  }
}
