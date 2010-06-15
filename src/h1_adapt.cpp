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
#include "element_to_refine.h"
#include "adapt.h"
#include "h1_adapt.h"

using namespace std;

H1Adapt::H1Adapt(const Tuple<Space*>& spaces) : Adapt(spaces) 
{
  for (int i = 0; i < spaces.size(); i++) {
    for (int j = 0; j < spaces.size(); j++) {
      if (i == j) {
        form[i][j] = h1_form<double, scalar>;
        ord[i][j]  = h1_form<Ord, Ord>;
      }
    }
  }
}

H1Adapt::H1Adapt(Space* s) : Adapt(Tuple<Space*>(s))
{
  int i = 0, j = 0;
  form[i][j] = h1_form<double, scalar>;
  ord[i][j]  = h1_form<Ord, Ord>;
} 

