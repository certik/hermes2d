#ifndef __HERMES2D_ADAPT_RESULT_H
#define __HERMES2D_ADAPT_RESULT_H

/// \brief A result of adaptivity step.
struct PUBLIC_API AdaptResult {
  int element_id; ///< Element that was split or modified.
  int split_type; ///< Split type.
  
  int sons_cnt; ///< A number of sons. Can be one of {1, 2, 4}.
  short sons_orders[4]; ///< Sons IDs

  AdaptResult() {};
  AdaptResult(int element_id, int split_type, const int* sons_orders);
};

#endif
