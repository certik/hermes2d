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
#include "limit_order.h"
#include "forms.h"
#include "refmap.h"
#include "integrals_h1.h"

#include "shapeset_h1_all.h"
#include "shapeset_l2_all.h"
#include "shapeset_hc_all.h"

#include "element_to_refine.h"
#include "ref_selectors/optimum_selector.h"
#include "ref_selectors/proj_based_selector.h"
#include "adapt.h"

#include "ref_selectors/h1_proj_based_selector.h"
#include "h1_adapt.h"
#include "ref_selectors/l2_proj_based_selector.h"
#include "l2_adapt.h"
#include "ref_selectors/hcurl_proj_based_selector.h"
#include "hcurl_adapt.h"
#include "python_adapt.h"

using namespace std;
using namespace RefinementSelectors;

PythonAdapt::~PythonAdapt() {
  if (adapt_instance != NULL)
    delete adapt_instance;
}

CandList PythonAdapt::convert_old_to_candlist(int adapt_type, bool iso_only) {
  if (iso_only) {
    switch (adapt_type) {
      case 0: return H2D_HP_ISO;
      case 1: return H2D_H_ISO;
      case 2: return H2D_P_ISO;
    }
  }
  else {
    switch (adapt_type) {
      case 0: return H2D_HP_ANISO;
      case 1: return H2D_H_ANISO;
      case 2: return H2D_P_ANISO;
    }
  }
  error("Unknown adaptity type %d", adapt_type);
  return H2D_HP_ANISO;
}

void PythonAdapt::set_biform(int i, int j, Adapt::biform_val_t bi_form, Adapt::biform_ord_t bi_ord) {
  assert_msg(adapt_instance != NULL, "Internal error: Instance is NULL.");
  adapt_instance->set_biform(i, j, bi_form, bi_ord);
}

double PythonAdapt::calc_error(MeshFunction* sln, MeshFunction* rsln) {
  return calc_error_n(1, sln, rsln);
}

double PythonAdapt::calc_error_2(MeshFunction* sln1, MeshFunction* sln2, MeshFunction* rsln1, MeshFunction* rsln2) {
  return calc_error_n(2, sln1, sln2, rsln1, rsln2);
}

double PythonAdapt::calc_error_n(int num, ...) {
  assert_msg(adapt_instance != NULL, "Internal error: Instance is NULL.");

  //prepare space for values
  Tuple<Solution*> sln, rsln;
  sln.reserve(num);
  rsln.reserve(num);

  //obtain values
  va_list ap;
  va_start(ap, num);
  for (int i = 0; i < num; i++)
    sln.push_back(va_arg(ap, Solution*)); //unsafe C-style type-casting to a child class: required by Python wrappers.
  for (int i = 0; i < num; i++)
    rsln.push_back(va_arg(ap, Solution*)); //unsafe C-style type-casting to a child class: required by Python wrappers.
  va_end(ap);

  //call error evaluation
  adapt_instance->set_solutions(sln, rsln);
  double error = adapt_instance->calc_error(H2D_TOTAL_ERROR_REL | H2D_ELEMENT_ERROR_ABS);
  return error;
}

H1OrthoHP::H1OrthoHP(int num, ...) {
  //obtain list of spaces
  Tuple<Space*> spaces;
  va_list ap;
  va_start(ap, num);
  for(int i = 0; i < num; i++)
    spaces.push_back(va_arg(ap, Space*));
  va_end(ap);

  //create the instance
  adapt_instance = new H1Adapt(spaces);
}

bool H1OrthoHP::adapt(double thr, int strat, int adapt_type, bool iso_only, int regularize, double conv_exp, int max_order, bool same_orders, double to_be_processed) {
  assert_msg(adapt_instance != NULL, "Internal error: Instance is NULL.");

  //create selector
  H1ProjBasedSelector selector(convert_old_to_candlist(adapt_type, iso_only), conv_exp, max_order);

  //call adapt
  bool result = adapt_instance->adapt(&selector, thr, strat, regularize, same_orders, to_be_processed);

  return result;
}

L2OrthoHP::L2OrthoHP(int num, ...) {
  //obtain list of spaces
  Tuple<Space*> spaces;
  va_list ap;
  va_start(ap, num);
  for(int i = 0; i < num; i++)
    spaces.push_back(va_arg(ap, Space*));
  va_end(ap);

  //create the instance
  adapt_instance = new L2Adapt(spaces);
}

bool L2OrthoHP::adapt(double thr, int strat, int adapt_type, bool iso_only, int regularize, double conv_exp, int max_order, bool same_orders, double to_be_processed) {
  assert_msg(adapt_instance != NULL, "Internal error: Instance is NULL.");

  //create selector
  L2ProjBasedSelector selector(convert_old_to_candlist(adapt_type, iso_only), conv_exp, max_order);

  //call adapt
  bool result = adapt_instance->adapt(&selector, thr, strat, regularize, same_orders, to_be_processed);

  return result;
}
