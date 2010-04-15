#include "trans.h"
#define countof(a) (sizeof(a)/sizeof(a[0]))

// taken form RefMap.cpp
extern H1ShapesetBeuchler ref_map_shapeset;

// vertical edge
double2 e0_pt[] = {
	{ -1,  1 },
	{ -1, -1 }
};

// horz edge
double2 e1_pt[] = {
	{  1,   -1 },
	{ -1,   -1 }
};

// last edge
double2 e2_pt[] = {
	{   -1,  1 },
	{ -0.5,  0.5 },
	{    0,  0 },
	{  0.5, -0.5 },
	{    1, -1 }
};

// modified get_phys_x,y
// @param[in] pt - point to tranforms
// @return transformed points
double2 *transform_element(Element *e, int np, double2 *pt)
{
	Quad2D *quad_2d = &g_quad_2d_std;

	ref_map_shapeset.set_mode(e->get_mode());

	// transform all x coordinates of the integration points
	double2 *tpt = new double2[np];
	memset(tpt, 0, np * sizeof(double2));

	int indices[70];

	// taken from refmap::set_active_element
	// prepare the shapes and coefficients of the reference map
	int j, k = 0;
	for (unsigned int i = 0; i < e->nvert; i++)
		indices[k++] = ref_map_shapeset.get_vertex_index(i);

	int o = e->cm->order;
	for (unsigned int i = 0; i < e->nvert; i++)
		for (j = 2; j <= o; j++)
			indices[k++] = ref_map_shapeset.get_edge_index(i, 0, j);

	if (e->is_quad()) o = make_quad_order(o, o);
	memcpy(indices + k, ref_map_shapeset.get_bubble_indices(o),
		ref_map_shapeset.get_num_bubbles(o) * sizeof(int));

	// taken from RefMap::get_phys_x
	for (int i = 0; i < e->cm->nc; i++)
	{
		for (j = 0; j < np; j++) {
			double fn = ref_map_shapeset.get_fn_value(indices[i], pt[j][0], pt[j][1], 0);
			tpt[j][0] += e->cm->coefs[i][0] * fn;
			tpt[j][1] += e->cm->coefs[i][1] * fn;
		}
	}

	return tpt;
}

double2 *transform(Element *e)
{
	double2 *tpt1 = NULL;
	double2 *tpt2 = NULL;
	double2 *tpt3 = NULL;

/*
	tpt = transform_element(e, countof(e0_pt), e0_pt);
	printf("edge #0\n");
	for (int i = 0; i < countof(e0_pt); i++)
		printf("%lf, %lf\n", tpt[i][0], tpt[i][1]);

	tpt = transform_element(e, countof(e1_pt), e1_pt);
	printf("edge #1\n");
	for (int i = 0; i < countof(e1_pt); i++)
		printf("%lf, %lf\n", tpt[i][0], tpt[i][1]);
*/

	tpt3 = transform_element(e, countof(e2_pt), e2_pt);
	//printf("edge #2\n");
	//for (int i = 0; i < countof(e2_pt); i++)
	//	printf("%lf, %lf\n", tpt3[i][0], tpt3[i][1]);
	return tpt3;
}

void element_polygonal_boundary(Element *e, double2 **tp, int *n)
{
	double2 *pt;

    //*tp = transform_element(e, countof(e2_pt), e2_pt);
    *n = 4;
    *tp = new double2[*n];
    pt = *tp;
    pt[0][0] = 1.;
    pt[0][1] = 1.;
    pt[1][0] = 2.;
    pt[1][1] = 2.;
    pt[2][0] = 3.;
    pt[2][1] = 3.;
    pt[3][0] = 4.;
    pt[3][1] = 4.;
}
