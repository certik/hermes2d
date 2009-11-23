#include <stdio.h>
#include "hermes2d.h"
#include "adapt_result.h"

AdaptResult::AdaptResult(int element_id, int split_type, const int* sons_orders)
  : element_id(element_id), split_type(split_type) {
  //get sons cnt
  if (split_type < 0)
    sons_cnt = 1;
  else if (split_type == 0)
    sons_cnt = 4;
  else if (split_type == 1 || split_type == 2)
    sons_cnt = 2;
  else {
    debug_log("E invalid or unsupported split type (%d) (AdaptResult::AdaptResult)\n", split_type);
    assert(false);
  }

  //copy orders
  for(int i = 0; i < sons_cnt; i++)
    this->sons_orders[i] = (short)sons_orders[i];
};
