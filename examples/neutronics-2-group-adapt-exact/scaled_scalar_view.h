#include "hermes2d.h"

class ScaledScalarView : public ScalarView 
{
	public:
		ScaledScalarView(const char* title = "ScaledScalarView", DEFAULT_WINDOW_POS) : 
			ScalarView(title, x, y, width, height) { }
		
		void scale(double sc = 1e3) {
			lin.lock_data();
			if (mode3d)
        yscale *= sc;
      else if (contours)
	      cont_step *= sc;
		  lin.unlock_data();   		  
		  refresh();
		}
};
