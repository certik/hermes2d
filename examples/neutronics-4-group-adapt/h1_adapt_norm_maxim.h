#ifndef __H1_ADAPT_NORM_MAXIM
#define __H1_ADAPT_NORM_MAXIM

enum NormType {
  NORM_EUCLEDIAN,
  NORM_MAXIMA
};

class H1AdaptNormMaxim : public H1Adapt {
  const NormType norm_type;

protected:
  double norm_fn_inf(MeshFunction* sln, RefMap* ru); ///< Function used to calculate L_oo norm of the solution.
  virtual scalar eval_norm(biform_val_t bi_fn, biform_ord_t bi_ord, MeshFunction *rsln1, MeshFunction *rsln2, RefMap *rrv1, RefMap *rrv2);

public:
  H1AdaptNormMaxim(NormType norm_type, const Tuple<Space*>& spaces) : H1Adapt(spaces), norm_type(norm_type) {};

};

#endif
