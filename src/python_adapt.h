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

#ifndef __H2D_PYTHON_ADAPT_H
#define __H2D_PYTHON_ADAPT_H

/// Python wrapper support class.
/** This class is not intended for use in C/C++ application. It is provided
 *  just for backward compatibility with Python wrappers. A regular C/C++ application
 *  should use children of the class Adapt instead. */
class H2D_API PythonAdapt {
protected:
  Adapt* adapt_instance; ///< Internal instance of adapt.

  /// Converts the old version of candidate generation to the new one.
  /** \param[in] adapt_type A type of adaptivity (0 = HP, 1 = H, 2 = P).
   *  \param[in] iso_only True if iso-only candidates has to be generated.
   *  \return A value identifying a predefined candidate types. */
  static RefinementSelectors::CandList convert_old_to_candlist(int adapt_type, bool iso_only);
public:
  PythonAdapt() : adapt_instance(NULL) {}; ///< Constructor.
  virtual ~PythonAdapt(); ///< Destructor.

  /// Sets user defined bilinear form which is used to calculate error. This method is just a wrapper. For details see Adapt::set_biform().
  void set_biform(int i, int j, Adapt::biform_val_t bi_form, Adapt::biform_ord_t bi_ord);

  double calc_error(MeshFunction* sln, MeshFunction* rsln); ///< A wrapper to Adapt::calc_error().
  double calc_error_2(MeshFunction* sln1, MeshFunction* sln2, MeshFunction* rsln1, MeshFunction* rsln2); ///< A wrapper to Adapt::calc_error().
  double calc_error_n(int n, ...); ///< A wrapper to Adapt::calc_error().
};

/// Python wrapper support class for H1 adaptivity.
/** This class is not intended for use in C/C++ application. It is provided
 *  just for backward compatibility with Python wrappers. A regular C/C++ application
 *  should use the class H1Adapt instead. */
class H2D_API H1OrthoHP : public PythonAdapt {
public:
  H1OrthoHP(int num, ...); ///< Contructor.

  /// A wrapper for the method Adapt::adapt().
  bool adapt(double thr, int strat = 0, int adapt_type = 0, bool iso_only = false, int regularize = -1,
             double conv_exp = 1.0, int max_order = -1, bool same_orders = false, double to_be_processed = 0.0);

};

/// Python wrapper support class for L2 adaptivity.
/** This class is not intended for use in C/C++ application. It is provided
 *  just for backward compatibility with Python wrappers. A regular C/C++ application
 *  should use the class L2Adapt instead. */
class H2D_API L2OrthoHP : public PythonAdapt {
public:
  L2OrthoHP(int num, ...); ///< Contructor.

  /// A wrapper for the method Adapt::adapt().
  bool adapt(double thr, int strat = 0, int adapt_type = 0, bool iso_only = false, int regularize = -1,
             double conv_exp = 1.0, int max_order = -1, bool same_orders = false, double to_be_processed = 0.0);

};

#endif
