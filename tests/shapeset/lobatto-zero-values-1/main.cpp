#include "hermes2d.h"
#define ERROR_SUCCESS                               0
#define ERROR_FAILURE                               -1

int P_INIT = 10;
double EPS = 10e-14;
// This test testing zero lobatto shape functions 
// values at vertices, on edges and bubbles.

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  H2DReader mloader;
  // We load the mesh on a (-1, 1)^2 domain.
  mloader.load("ref_square.mesh", &mesh);            

  // Create an H1 space with default shapeset,
  // natural BC, and linear elements.
  H1Space space(&mesh, NULL, NULL, P_INIT);
  // The type of element, mesh_mode = 4 means a rectangle element.
  int mesh_mode = 4;
  int n = get_num_dofs(&space);

  int *fn_idx = new int [n];
  int m = 0;
  int order = P_INIT;
  int vertex1 = 0;
  int vertex2 = 1;
  int vertex3 = 2;
  int vertex4 = 3;

  double x = 0.0, y = 0.0, value = 0.0;
  double x1 = -1.0, y1 = -1.0;
  double x2 =  1.0, y2 = -1.0;
  double x3 =  1.0, y3 =  1.0;
  double x4 = -1.0, y4 =  1.0;

  info("Testing................");
  // Check vertex functions.
  info("Check vertex functions.");
  for (int i = 0; i < mesh_mode; i++, m++)
  {
    fn_idx[m] = space.get_shapeset()->get_vertex_index(i);
    info("Check vertex function [%d]", m);

    if (i == vertex1)
    {
      // Vertices.
      value = space.get_shapeset()->get_fn_value(fn_idx[i], x2, y2, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }

      value = space.get_shapeset()->get_fn_value(fn_idx[i], x3, y3, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
      value = space.get_shapeset()->get_fn_value(fn_idx[i], x4, y4, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
      printf("Vertices... ok\n");

      // Edges.
      printf("Edge  2.");
      for (int j = 0; j < order-1; j++)
      {
        x = x3 - (j+1)*(x3 - x2)/order;
        y = y3 - (j+1)*(y3 - y2)/order;
        value = space.get_shapeset()->get_fn_value(fn_idx[i], x, y, 0);
        if (value >= EPS)
        {
          printf("Failed\n");
          return ERROR_FAILURE;
        }
      }
      printf("... ok\n");
      printf("Edge  3.");
      for (int j = 0; j < order-1; j++)
      {
        x = x3 - (j+1)*(x3 - x4)/order;
        y = y3 - (j+1)*(y3 - y4)/order;
        value = space.get_shapeset()->get_fn_value(fn_idx[i], x, y, 0);
        if (value >= EPS)
        {
          printf("Failed\n");
          return ERROR_FAILURE;
        }
      }
      printf("... ok\n");
    }    

    if (i == vertex2)
    {
      // Vertices.
      value = space.get_shapeset()->get_fn_value(fn_idx[i], x1, y1, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
      value = space.get_shapeset()->get_fn_value(fn_idx[i], x3, y3, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
      value = space.get_shapeset()->get_fn_value(fn_idx[i], x4, y4, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
      printf("Vertices... ok\n");

      // Edges.
      printf("Edge  4.");
      for (int j = 0; j < order-1; j++)
      {
        x = x4 - (j+1)*(x4 - x1)/order;
        y = y4 - (j+1)*(y4 - y1)/order;
        value = space.get_shapeset()->get_fn_value(fn_idx[i], x, y, 0);
        if (value >= EPS)
        {
          printf("Failed\n");
          return ERROR_FAILURE;
        }
      }
      printf("... ok\n");
      printf("Edge  3.");
      for (int j = 0; j < order-1; j++)
      {
        x = x4 - (j+1)*(x4 - x3)/order;
        y = y4 - (j+1)*(y4 - y3)/order;
        value = space.get_shapeset()->get_fn_value(fn_idx[i], x, y, 0);
        if (value >= EPS)
        {
          printf("Failed\n");
          return ERROR_FAILURE;
        }
      }
      printf("... ok\n");
    }    

    if (i == vertex3)
    {
      // Vertices.
      value = space.get_shapeset()->get_fn_value(fn_idx[i], x1, y1, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
      value = space.get_shapeset()->get_fn_value(fn_idx[i], x2, y2, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
      value = space.get_shapeset()->get_fn_value(fn_idx[i], x4, y4, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
      printf("Vertices... ok\n");

      // Edges.
      printf("Edge  4.");
      for (int j = 0; j < order-1; j++)
      {
        x = x1 - (j+1)*(x1 - x4)/order;
        y = y1 - (j+1)*(y1 - y4)/order;
        value = space.get_shapeset()->get_fn_value(fn_idx[i], x, y, 0);
        if (value >= EPS)
        {
          printf("Failed\n");
          return ERROR_FAILURE;
        }
      }
      printf("... ok\n");
      printf("Edge  1.");
      for (int j = 0; j < order-1; j++)
      {
        x = x1 - (j+1)*(x1 - x2)/order;
        y = y1 - (j+1)*(y1 - y2)/order;
        value = space.get_shapeset()->get_fn_value(fn_idx[i], x, y, 0);
        if (value >= EPS)
        {
          printf("Failed\n");
          return ERROR_FAILURE;
        }
      }
      printf("... ok\n");
    }    

    if (i == vertex4)
    {
      // Vertices.
      value = space.get_shapeset()->get_fn_value(fn_idx[i], x1, y1, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
      value = space.get_shapeset()->get_fn_value(fn_idx[i], x2, y2, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
      value = space.get_shapeset()->get_fn_value(fn_idx[i], x3, y3, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
      printf("Vertices... ok\n");

      // Edges.
      printf("Edge  2.");
      for (int j = 0; j < order-1; j++)
      {
        x = x2 - (j+1)*(x2 - x3)/order;
        y = y2 - (j+1)*(y2 - y3)/order;
        value = space.get_shapeset()->get_fn_value(fn_idx[i], x, y, 0);
        if (value >= EPS)
        {
          printf("Failed\n");
          return ERROR_FAILURE;
        }
      }
      printf("... ok\n");
      printf("Edge  1.");
      for (int j = 0; j < order-1; j++)
      {
        x = x2 - (j+1)*(x2 - x1)/order;
        y = y2 - (j+1)*(y2 - y1)/order;
        value = space.get_shapeset()->get_fn_value(fn_idx[i], x, y, 0);
        if (value >= EPS)
        {
          printf("Failed\n");
          return ERROR_FAILURE;
        }
      }
      printf("... ok\n");
    }   
  }

  // Check edge functions.
  info("Check edge functions.");
  for (int edge_order = 2; edge_order <= order; edge_order++)
  {
    for (int j = 0; j < mesh_mode; j++, m++)
    {
      fn_idx[m] = space.get_shapeset()->get_edge_index(j, 0, edge_order);
      info("Check edge function [%d]", m);

      // Vertices.
      value = space.get_shapeset()->get_fn_value(fn_idx[m], x1, y1, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
      value = space.get_shapeset()->get_fn_value(fn_idx[m], x2, y2, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
      value = space.get_shapeset()->get_fn_value(fn_idx[m], x3, y3, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
      value = space.get_shapeset()->get_fn_value(fn_idx[m], x4, y4, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
      printf("Vertices... ok\n");

      // Edge 1.
      if (m%4 == 0)
      {
        printf("Egde  1.");
        // Check edge 2.
        for (int j = 0; j < order-1; j++)
        {
          x = x2;
          y = y2 - (j+1)*(y2 - y3)/order;
          value = space.get_shapeset()->get_fn_value(fn_idx[m], x, y, 0);
          if (value >= EPS)
          {
            printf("Failed\n");
            return ERROR_FAILURE;
          }
        }
        // Check edge 3.
        for (int j = 0; j < order-1; j++)
        {
          x = x3;
          y = y3 - (j+1)*(y3 - y4)/order;
          value = space.get_shapeset()->get_fn_value(fn_idx[m], x, y, 0);
          if (value >= EPS)
          {
            printf("Failed\n");
            return ERROR_FAILURE;
          }
        }
        // Check edge 4.
        for (int j = 0; j < order-1; j++)
        {
          x = x1 - (j+1)*(x1 - x4)/order;
          y = y1;
          value = space.get_shapeset()->get_fn_value(fn_idx[m], x, y, 0);
          if (value >= EPS)
          {
            printf("Failed\n");
            return ERROR_FAILURE;
          }
        }
        printf("... ok\n");
      }

      // Edge 2.
      if (m%4 == 1)
      {
        printf("Egde  2.");
        // Check edge 1.
        for (int j = 0; j < order-1; j++)
        {
          x = x2 - (j+1)*(x2 - x1)/order;
          y = y2;
          value = space.get_shapeset()->get_fn_value(fn_idx[m], x, y, 0);
          if (value >= EPS)
          {
            printf("Failed\n");
            return ERROR_FAILURE;
          }
        }
        // Check edge 3.
        for (int j = 0; j < order-1; j++)
        {
          x = x3 - (j+1)*(x3 - x4)/order;;
          y = y3;
          value = space.get_shapeset()->get_fn_value(fn_idx[m], x, y, 0);
          if (value >= EPS)
          {
            printf("Failed\n");
            return ERROR_FAILURE;
          }
        }
        // Check edge 4.
        for (int j = 0; j < order-1; j++)
        {
          x = x4;
          y = y4 - (j+1)*(x4 - x1)/order;
          value = space.get_shapeset()->get_fn_value(fn_idx[m], x, y, 0);
          if (value >= EPS)
          {
            printf("Failed\n");
            return ERROR_FAILURE;
          }
        }
        printf("... ok\n");
      }

      // Edge 3.
      if (m%4 == 2)
      {
        printf("Egde  3.");
        // Check edge 2.
        for (int j = 0; j < order-1; j++)
        {
          x = x2;
          y = y2 - (j+1)*(y2 - y3)/order;
          value = space.get_shapeset()->get_fn_value(fn_idx[m], x, y, 0);
          if (value >= EPS)
          {
            printf("Failed\n");
            return ERROR_FAILURE;
          }
        }
        // Check edge 4.
        for (int j = 0; j < order-1; j++)
        {
          x = x1;
          y = y1 - (j+1)*(y1 - y4)/order;
          value = space.get_shapeset()->get_fn_value(fn_idx[m], x, y, 0);
          if (value >= EPS)
          {
            printf("Failed\n");
            return ERROR_FAILURE;
          }
        }
        // Check edge 1.
        for (int j = 0; j < order-1; j++)
        {
          x = x1 - (j+1)*(x1 - x2)/order;
          y = y1;
          value = space.get_shapeset()->get_fn_value(fn_idx[m], x, y, 0);
          if (value >= EPS)
          {
            printf("Failed\n");
            return ERROR_FAILURE;
          }
        }
        printf("... ok\n");
      }

      // Edge 4.
      if (m%4 == 3)
      {
        printf("Egde  4.");
        // Check edge 1.
        for (int j = 0; j < order-1; j++)
        {
          x = x2 - (j+1)*(x2 - x1)/order;
          y = y2;
          value = space.get_shapeset()->get_fn_value(fn_idx[m], x, y, 0);
          if (value >= EPS)
          {
            printf("Failed\n");
            return ERROR_FAILURE;
          }
        }
        // Check edge 2.
        for (int j = 0; j < order-1; j++)
        {
          x = x2;
          y = y2 - (j+1)*(y2 - y3)/order;
          value = space.get_shapeset()->get_fn_value(fn_idx[m], x, y, 0);
          if (value >= EPS)
          {
            printf("Failed\n");
            return ERROR_FAILURE;
          }
        }
        // Check edge 3.
        for (int j = 0; j < order-1; j++)
        {
          x = x3 - (j+1)*(x3 - x4)/order;
          y = y3;
          value = space.get_shapeset()->get_fn_value(fn_idx[m], x, y, 0);
          if (value >= EPS)
          {
            printf("Failed\n");
            return ERROR_FAILURE;
          }
        }
        printf("... ok\n");
      }
    }
  }

  // Check bubble functions.
  info("Check bubble functions.");
  int number_bubble = space.get_shapeset()->get_num_bubbles(H2D_MAKE_QUAD_ORDER(order, order));
  int *bubble_idx = space.get_shapeset()->get_bubble_indices(H2D_MAKE_QUAD_ORDER(order, order));
  for (int i = 0; i < number_bubble; i++, m++ )
  {
    fn_idx[m] = bubble_idx[i];
    info("Check bubble function [%d]", m);

    printf("Egde 1 and vertex 1.");
    // Check edge 1 and vertex 1.
    for (int j = 0; j < order; j++)
    {
      x = x1 - (j)*(x2 - x1)/order;
      y = y1;
      value = space.get_shapeset()->get_fn_value(fn_idx[m], x, y, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
    }
    printf("... ok\n");

    printf("Egde 2 and vertex 2.");
    // Check edge 2 and vertex 2.
    for (int j = 0; j < order; j++)
    {
      x = x2;
      y = y2 - (j)*(y2 - y3)/order;
      value = space.get_shapeset()->get_fn_value(fn_idx[m], x, y, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
    }
    printf("... ok\n");

    printf("Egde 3 and vertex 3.");
    // Check edge 3 and vertex 3.
    for (int j = 0; j < order; j++)
    {
      x = x3 - (j)*(x3 - x4)/order;
      y = y3;
      value = space.get_shapeset()->get_fn_value(fn_idx[m], x, y, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
    }
    printf("... ok\n");

    printf("Egde 4 and vertex 4.");
    // Check edge 4 and vertex 4.
    for (int j = 0; j < order; j++)
    {
      x = x4;
      y = y4 - (j)*(y4 - y4)/order;
      value = space.get_shapeset()->get_fn_value(fn_idx[m], x, y, 0);
      if (value >= EPS)
      {
        printf("Failed\n");
        return ERROR_FAILURE;
      }
    }
    printf("... ok\n");
  }
  info("ndof = %d", n);

  printf("Success!\n");
 
  delete [] fn_idx;
  return ERROR_SUCCESS;
}

