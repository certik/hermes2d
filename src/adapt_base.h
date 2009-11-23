#ifndef __HERMES2D_ADAPT_BASE_L2_H
#define __HERMES2D_ADAPT_BASE_L2_H

/// \brief A result of adaptivity step.
struct PUBLIC_API AdaptResult {
  int element_id; ///< Element that was split or modified.
  int sons_cnt; ///< A number of sons. Can be one of {1, 2, 4}.
  int max_sons_order; ///< A maximum order of sons.
  AdaptResult() {};
  AdaptResult(int element_id, int sons_cnt, int max_sons_order) : element_id(element_id), sons_cnt(sons_cnt), max_sons_order(max_sons_order) {};
};

#endif
