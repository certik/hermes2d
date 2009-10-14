/// THIS CLASS IS DEPRECATED

#include "common.h"
#include "discrete.h"
#include "space.h"
#include "precalc.h"
#include "refmap.h"
#include <umfpack.h>
#include "quad_all.h"
#include "traverse.h"


//// UMFPACK and other stuff ///////////////////////////////////////////////////////////////////////

#ifndef COMPLEX
  // real case: no changes in calling UMFPACK
  #define umfpack_symbolic(m, n, Ap, Ai, Ax, S, C, I)   umfpack_di_symbolic(m, n, Ap, Ai, Ax, S, C, I)
  #define umfpack_numeric(Ap, Ai, Ax, S, N, C, I)       umfpack_di_numeric(Ap, Ai, Ax, S, N, C, I)
  #define umfpack_solve(sys, Ap, Ai, Ax, X, B, N, C, I) umfpack_di_solve(sys, Ap, Ai, Ax, X, B, N, C, I)
  #define umfpack_free_symbolic                         umfpack_di_free_symbolic
  #define umfpack_free_numeric                          umfpack_di_free_numeric
  #define umfpack_defaults                              umfpack_di_defaults
#else
  // macros for calling complex UMFPACK in packed-complex mode
  #define umfpack_symbolic(m, n, Ap, Ai, Ax, S, C, I) \
          umfpack_zi_symbolic(m, n, Ap, Ai, (double*) (Ax), NULL, S, C, I)
  #define umfpack_numeric(Ap, Ai, Ax, S, N, C, I) \
          umfpack_zi_numeric(Ap, Ai, (double*) (Ax), NULL, S, N, C, I)
  #define umfpack_solve(sys, Ap, Ai, Ax, X, B, N, C, I) \
          umfpack_zi_solve(sys, Ap, Ai, (double*) (Ax), NULL, (double*) (X), NULL, (double*) (B), NULL, N, C, I)
  #define umfpack_free_symbolic umfpack_zi_free_symbolic
  #define umfpack_free_numeric  umfpack_zi_free_numeric
  #define umfpack_defaults      umfpack_zi_defaults
#endif

static void umfpack_status(int status)
{
  switch (status)
  {
    case UMFPACK_OK:                            info ("UMFPACK: OK status"); break;
    case UMFPACK_WARNING_singular_matrix:
      warn ("UMFPACK: singular stiffness matrix!"); break;
    case UMFPACK_ERROR_out_of_memory:           error("UMFPACK: out of memory!");
    case UMFPACK_ERROR_argument_missing:        error("UMFPACK: argument missing");
    case UMFPACK_ERROR_invalid_Symbolic_object: error("UMFPACK: invalid Symbolic object");
    case UMFPACK_ERROR_invalid_Numeric_object:  error("UMFPACK: invalid Numeric object");
    case UMFPACK_ERROR_different_pattern:       error("UMFPACK: different pattern");
    case UMFPACK_ERROR_invalid_system:          error("UMFPACK: invalid system");
    case UMFPACK_ERROR_n_nonpositive:           error("UMFPACK: n nonpositive");
    case UMFPACK_ERROR_invalid_matrix:          error("UMFPACK: invalid matrix");
    case UMFPACK_ERROR_internal_error:          error("UMFPACK: internal error");
    default:                                    error("UMFPACK: unknown error");
  }
}


/*static int default_order_table_tri[] =
{
  0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16,
  17, 18, 19, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
  20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
  20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
  20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20
};

static int default_order_table_quad[] =
{
  1, 1, 3, 3, 5, 5, 7, 7, 9, 9, 11, 11, 13, 13, 15, 15, 17,
  17, 19, 19, 21, 21, 23, 23, 24, 24, 24, 24, 24, 24, 24,
  24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
  24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24,
  24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24
};

int  g_max_order;
int* g_order_table_quad = default_order_table_quad;
int* g_order_table_tri  = default_order_table_tri;
int* g_order_table = NULL;
bool warned_order = false;*/
extern bool warned_order;

void qsort_int(int* pbase, size_t total_elems); // defined in qsort.cpp


//// precalculation of nonzero element positions //////////////////////////////////////////////////

// How it works: a special assembly-like procedure is invoked before the real assembly, whose goal is
// to determine the positions of nonzero elements in the stiffness matrix. Naturally, the bilinear
// form is not being evaluated at this point, just the global DOF indices are used (the array 'dof').
// The nonzero positions are simply accumulated for each row(*) in an array-like structure. Because
// the lengths of the arrays are not known ahead, they are allocated in blocks called pages (see the
// structure Page below). Any time the array is full and a new element position (index)
// must be added, a new page is allocated. Because of the nature of the element-by-element assembly,
// it is probable that one index will be inserted more than once, which corresponds to adding several
// values to a single matrix entry. The array thus has to be sorted at the end, which allows counting
// of the nonzero positions while disregarding their duplicities. The duplicity for each matrix position
// is about two on average, hence the precalculation process requires about two thirds of the total
// memory that will be required for the final matrix (8 bytes vs. 12 bytes per nonzero element).
// After counting the nonzero elements all pages are freed, so no matrix memory is wasted.
//
// (*) For UMFPACK column positions are accumulated.

static const int page_size = 62;

struct Page
{
  int count;
  int idx[page_size];
  Page* next;
};


static inline void page_add_ij(Page** pages, int i, int j)
{
  if (pages[i] == NULL || pages[i]->count >= page_size)
  {
    Page* new_page = new Page;
    new_page->count = 0;
    new_page->next = pages[i];
    pages[i] = new_page;
  }
  pages[i]->idx[pages[i]->count++] = j;
}


void DiscreteProblem::precalculate_sparse_structure(Page** pages)
{
  int i, j, m, n;
  AsmList* al = new AsmList[neq];

  // init multi-mesh traversal
  Mesh** meshes = new Mesh*[neq];
  for (i = 0; i < neq; i++)
    meshes[i] = spaces[i]->get_mesh();
  Traverse trav;
  trav.begin(neq, meshes);

  // loop through all triangles
  Element** e;
  while ((e = trav.get_next_state(NULL, NULL)) != NULL)
  {
    // obtain assembly lists for the element at all spaces
    for (i = 0; i < neq; i++)
      spaces[i]->get_element_assembly_list(e[i], al + i);
      // todo: neziskavat znova, pokud se element nezmenil

    // go through all equation-blocks of the local stiffness matrix
    for (m = 0; m < neq; m++)
    {
      for (n = 0; n < neq; n++)
      {
        BiForm* bf = biform[m] + n;
        if (bf->sym == NULL && bf->unsym == NULL && bf->surf == NULL) continue;

        // pretend assembling of the element stiffness matrix
        for (j = 0; j < al[n].cnt; j++)
        {
          // skip dirichlet dofs in 'j'
          if (al[n].dof[j] < 0) continue;

          for (i = 0; i < al[m].cnt; i++)
          {
            // skip dirichlet dofs in 'i'
            if (al[m].dof[i] < 0) continue;

            // register the corresponding nonzero matrix element
            page_add_ij(pages, al[n].dof[j], al[m].dof[i]); // column-oriented (UMFPACK)
            //page_add_ij(pages, al[m].dof[i], al[n].dof[j]); // row-oriented (PETSc)
          }
        }
      }
    }
  }

  trav.finish();
  delete [] meshes;
  delete [] al;
}


static int sort_and_store_indices(Page* page, int* buffer, int* max)
{
  // gather all pages in the buffer, deleting them along the way
  int* end = buffer;
  while (page != NULL)
  {
    memcpy(end, page->idx, sizeof(int) * page->count);
    end += page->count;
    Page* tmp = page;
    page = page->next;
    delete tmp;
  }

  // sort the indices and remove duplicities
  qsort_int(buffer, end - buffer);
  int *q = buffer;
  for (int *p = buffer, last = -1;  p < end;  p++)
    if (*p != last)
      *q++ = last = *p;

  return q - buffer;
}


static int get_num_indices(Page** pages, int ndofs)
{
  int total = 0;
  for (int i = 0; i < ndofs; i++)
    for (Page* page = pages[i]; page != NULL; page = page->next)
      total += page->count;

  return total;
}


//// matrix creation ///////////////////////////////////////////////////////////////////////////////

void DiscreteProblem::create_matrix()
{
  // remove any previous matrix
  free_matrix_indices();
  free_matrix_values();

  // calculate the total number of DOFs
  ndofs = 0;
  for (int i = 0; i < neq; i++)
    ndofs += spaces[i]->get_num_dofs();
  if (!quiet) verbose("Ndofs: %d", ndofs);
  if (!ndofs) return;

  // get row and column indices of nonzero matrix elements
  Page** pages = new Page*[ndofs];
  memset(pages, 0, sizeof(Page*) * ndofs);
  if (!quiet) { verbose("Calculating matrix sparse structure..."); begin_time(); }
  precalculate_sparse_structure(pages);

  // initialize the arrays Ap and Ai
  Ap = (int*) malloc(sizeof(int) * (ndofs+1));
  int aisize = get_num_indices(pages, ndofs);
  Ai = (int*) malloc(sizeof(int) * aisize);
  if (Ai == NULL) error("Out of memory. Could not allocate the array Ai.");

  // sort the indices and remove duplicities, insert into Ai
  int i, pos = 0, num;
  for (i = 0; i < ndofs; i++)
  {
    Ap[i] = pos;
    pos += sort_and_store_indices(pages[i], Ai + pos, Ai + aisize);
  }
  Ap[i] = pos;
  if (!quiet) verbose("  Nonzeros: %d\n  Total matrix size: %0.1lf MB\n  (time: %g sec)",
                         pos, (double) get_matrix_size() / (1024*1024), end_time());
  delete [] pages;

  // shrink Ai to the actual size
  int* oldAi = Ai;
  Ai = (int*) realloc(Ai, sizeof(int) * pos);
  if (oldAi != Ai) warn("Realloc moved Ai when shrinking."); // this should not happen

  // UMFPACK: perform symbolic analysis of the matrix
  if (!quiet) { verbose("Performing UMFPACK symbolic analysis..."); begin_time(); }
  int status = umfpack_symbolic(ndofs, ndofs, Ap, Ai, NULL, &Symbolic, NULL, NULL);
  if (status != UMFPACK_OK) umfpack_status(status);
  if (!quiet) verbose("  (time: %g sec)", end_time());

  equi = (double*) malloc(sizeof(double) * ndofs);
  if (equi == NULL) error("Out of memory. Error allocating the equilibration vector.");
  for (int i = 0; i < ndofs; i++)
    equi[i] = 1.0;
  is_equi = false;
}


//// assembly //////////////////////////////////////////////////////////////////////////////////////

void DiscreteProblem::insert_matrix(scalar** mat, int* iidx, int* jidx, int ilen, int jlen)
{
  for (int j = 0; j < jlen; j++)
  {
    int col = jidx[j];
    if (col < 0) continue;
    int* cidx = Ai + Ap[col];
    int  clen = Ap[col+1] - Ap[col];
    scalar* cval = Ax + Ap[col];

    for (register int i = 0; i < ilen; i++)
    {
      register int row = iidx[i];
      if (row < 0) continue;
      register int lo = 0, hi = clen-1, mid;
      while (1)
      {
        mid = (lo + hi) >> 1;
        if (row < cidx[mid])
          hi = mid-1;
        else if (row > cidx[mid])
          lo = mid+1;
        else
          break;
        if (lo > hi)
          error("Out-dated sparse matrix structure! Spaces must not change between the calls to"
                " create_stiffness_matrix() and assemble_stiffness_matrix_and_rhs().");
      }
      //assert(cidx[mid] == idx[i]);
      cval[mid] += mat[i][j];
      //while (i < ilen-1 && mid < clen-1 && cidx[++mid] == iidx[++i]);
      //  cval[mid] += mat[i][j];
    }
  }
}


void DiscreteProblem::assemble_matrix_and_rhs(bool rhsonly)
{
  int i, j, k, l, m, n;
  bool bnd[4], nat[neq];
  EdgePos ep[4];

  if (!ndofs) return;
  warned_order = false;

  if (!rhsonly)
  {
    alloc_matrix_values();
    if (!quiet) { verbose("Assembling stiffness matrix..."); begin_time(); }
  }
  else
  {
    memset(RHS, 0, sizeof(scalar) * ndofs);
    if (!quiet) { verbose("Assembling RHS..."); begin_time(); }
  }

  // create slave pss's for test functions, init quadrature points
  PrecalcShapeset* spss[neq];
  PrecalcShapeset *fu, *fv;
  for (i = 0; i < neq; i++)
  {
    spss[i] = new PrecalcShapeset(pss[i]);
    pss [i]->set_quad_2d(&g_quad_2d_std);
    spss[i]->set_quad_2d(&g_quad_2d_std);
  }

  // initialize buffer
  buffer = NULL;
  mat_size = 0;
  get_matrix_buffer(9);

  // initialize assembly lists, refmap
  AsmList al[neq], *am, *an;

  RefMap* refmap = new RefMap[neq];
  for (i = 0; i < neq; i++)
    refmap[i].set_quad_2d(&g_quad_2d_std);
  for (i = 0; i < num_extern; i++)
    extern_fns[i]->set_quad_2d(&g_quad_2d_std);

  // init multi-mesh traversal
  int nm = neq + num_extern;
  Mesh* meshes[nm];
  Transformable* fn[nm];
  for (i = 0; i < neq; i++)
    meshes[i] = spaces[i]->get_mesh();
  memcpy(fn, pss, neq * sizeof(Transformable*));
  for (i = 0; i < num_extern; i++) {
    meshes[neq+i] = extern_fns[i]->get_mesh();
    fn[neq+i] = extern_fns[i];
  }
  // todo: kdyz maji nektere slozky stejnou sit, at sdili i refmapy
  //  - ale to bysme potrebovali slave RefMap

  // loop through all elements
  Element** e;
  Traverse trav;
  trav.begin(nm, meshes, fn);
  while ((e = trav.get_next_state(bnd, ep)) != NULL)
  {
    // set maximum integration order for use in integrals, see limit_order()
    update_limit_table(e[0]->get_mode());

    // obtain assembly lists for the element at all spaces, set appropriate mode for each pss
    for (i = 0; i < neq; i++)
    {
      spaces[i]->get_element_assembly_list(e[i], al + i);
      // todo: neziskavat znova, pokud se element nezmenil
      if (is_equi)
        for (j = 0; j < al[i].cnt; j++)
          if (al[i].dof[j] >= 0)
            al[i].coef[j] /= equi[al[i].dof[j]];

      spss[i]->set_active_element(e[i]);
      spss[i]->set_master_transform();
      refmap[i].set_active_element(e[i]);
      refmap[i].force_transform(pss[i]->get_transform(), pss[i]->get_ctm());
    }

    // go through all equation-blocks of the element stiffness matrix, assemble volume integrals
    for (m = 0, am = al; m < neq; m++, am++)
    {
      fv = spss[m];
      if (!rhsonly)
      {
        for (n = 0, an = al; n < neq; n++, an++)
        {
          fu = pss[n];
          BiForm* bf = biform[m] + n;
          if (!bf->sym && !bf->unsym) continue;
          if (bf->unsym == BF_SYM || bf->unsym == BF_ANTISYM) continue;
          bool tra = (biform[n][m].unsym == BF_SYM || biform[n][m].unsym == BF_ANTISYM);

          // assemble the (m,n)-block of the stiffness matrix
          scalar sy, un, **mat = get_matrix_buffer(std::max(am->cnt, an->cnt));
          for (i = 0; i < am->cnt; i++)
          {
            if (!tra && (k = am->dof[i]) < 0) continue;
            fv->set_active_shape(am->idx[i]);

            // unsymmetric block
            if (!bf->sym)
            {
              for (j = 0; j < an->cnt; j++)
              {
                fu->set_active_shape(an->idx[j]);
                un = bf->unsym(fu, fv, refmap+n, refmap+m) * an->coef[j] * am->coef[i];
                if (an->dof[j] < 0) Dir[k] -= un; else mat[i][j] = un;
              }
            }
            // symmetric block
            else
            {
              for (j = 0; j < an->cnt; j++)
              {
                scalar coef = an->coef[j] * am->coef[i];
                if (an->dof[j] < 0)
                {
                  fu->set_active_shape(an->idx[j]);
                  un = bf->unsym ? bf->unsym(fu, fv, refmap+n, refmap+m) * coef : 0.0;
                  sy = bf->sym(fu, fv, refmap+n, refmap+m) * coef;
                  Dir[k] -= (un + sy);
                }
                if (j >= i)
                {
                  fu->set_active_shape(an->idx[j]);
                  un = bf->unsym ? bf->unsym(fu, fv, refmap+n, refmap+m) * coef : 0.0;
                  mat[j][i] = sy = bf->sym  (fu, fv, refmap+n, refmap+m) * coef;
                  mat[i][j] = (un + sy);
                }
                else if (bf->unsym)
                {
                  fu->set_active_shape(an->idx[j]);
                  mat[i][j] += bf->unsym(fu, fv, refmap+n, refmap+m) * coef;
                }
              }
            }
          }

          // insert the local stiffness matrix into the global one
          insert_matrix(mat, am->dof, an->dof, am->cnt, an->cnt);

          // insert also the off-diagonal (anti-)symmetric block, if required
          if (tra)
          {
            if (biform[n][m].unsym == BF_ANTISYM) chsgn(mat, am->cnt, an->cnt);
            transpose(mat, am->cnt, an->cnt);
            insert_matrix(mat, an->dof, am->dof, an->cnt, am->cnt);

            // we also need to take care of the RHS...
            for (j = 0; j < am->cnt; j++)
              if (am->dof[j] < 0)
                for (i = 0; i < an->cnt; i++)
                  if (an->dof[i] >= 0)
                    Dir[an->dof[i]] -= mat[i][j];
          }
        }
      }

      // assemble rhs (linear form)
      if (!liform[m].lf) continue;
      for (i = 0; i < am->cnt; i++)
      {
        if (am->dof[i] < 0) continue;
        fv->set_active_shape(am->idx[i]);
        RHS[am->dof[i]] += liform[m].lf(fv, refmap+m) * am->coef[i];
      }
    }

    // assemble surface integrals now: loop through boundary edges of the element
    if (rhsonly) continue; // fixme
    for (int edge = 0; edge < e[0]->nvert; edge++)
    {
      if (!bnd[edge]) continue;

      // obtain the list of shape functions which are nonzero on this edge
      for (i = 0; i < neq; i++)
        if ((nat[i] = (spaces[i]->bc_type_callback(ep[edge].marker) == BC_NATURAL)))
          spaces[i]->get_edge_assembly_list(e[i], edge, al + i);

      // loop through the equation-blocks
      for (m = 0, am = al; m < neq; m++, am++)
      {
        if (!nat[m]) continue;
        fv = spss[m];
        ep[edge].base = trav.get_base();
        ep[edge].space_v = spaces[m];
        for (n = 0, an = al; n < neq; n++, an++)
        {
          if (!nat[n]) continue;
          BiForm* bf = biform[m] + n;
          if (!bf->surf) continue;
          fu = pss[n];
          ep[edge].space_u = spaces[n];

          // assemble the surface part of the bilinear form
          scalar bi, **mat = get_matrix_buffer(std::max(am->cnt, an->cnt));
          for (i = 0; i < am->cnt; i++)
          {
            if ((k = am->dof[i]) < 0) continue;
            fv->set_active_shape(am->idx[i]);
            for (j = 0; j < an->cnt; j++)
            {
              fu->set_active_shape(an->idx[j]);
              bi = bf->surf(fu, fv, refmap+n, refmap+m, ep+edge) * an->coef[j] * am->coef[i];
              if (an->dof[j] >= 0) mat[i][j] = bi; else Dir[k] -= bi;
            }
          }
          insert_matrix(mat, am->dof, an->dof, am->cnt, an->cnt);
        }

        // assemble the surface part of the linear form
        if (!liform[m].surf) continue;
        for (i = 0; i < am->cnt; i++)
        {
          if (am->dof[i] < 0) continue;
          fv->set_active_shape(am->idx[i]);
          RHS[am->dof[i]] += liform[m].surf(fv, refmap+m, ep+edge) * am->coef[i];
        }
      }
    }
  }

  trav.finish();
  for (i = 0; i < ndofs; i++)
    RHS[i] += Dir[i];

  if (!quiet) verbose("  (time: %g sec)", end_time());
  for (i = 0; i < neq; i++) delete spss[i];
  delete [] buffer;
  delete [] refmap;
}


int DiscreteProblem::get_matrix_size() const
{
  return (sizeof(int) + sizeof(scalar)) * Ap[ndofs] + sizeof(scalar) * 2 * ndofs;
}

void DiscreteProblem::alloc_matrix_values()
{
  free_matrix_values();

  Ax  = (scalar*) malloc(sizeof(scalar) * Ap[ndofs]);
  if (Ax == NULL) error("Out of memory. Error allocating stiffness matrix (Ax).");
  memset(Ax,  0, sizeof(scalar) * Ap[ndofs]);

  RHS = (scalar*) malloc(sizeof(scalar) * ndofs);
  Dir = (scalar*) malloc(sizeof(scalar) * ndofs);
  if (RHS == NULL || Dir == NULL) error("Out of memory. Error allocating the RHS vector.");
  memset(RHS, 0, sizeof(scalar) * ndofs);
  memset(Dir, 0, sizeof(scalar) * ndofs);
}


//// equilibration /////////////////////////////////////////////////////////////////////////////////

static inline int find_entry(int *Ai, int i, int lo, int hi)
{
  while(1)
  {
    int k = (lo + hi) >> 1;
    if (i < Ai[k])
      hi = k - 1;
    else if (i > Ai[k])
      lo = k + 1;
    else
      return k;
    if (lo > hi)
      return -1;
  }
}


void DiscreteProblem::equilibrate_matrix()
{
  int i, j;
  verbose("Equilibrating stiffness matrix...");

  // obtain the square root of the diagonal
  for (i = 0; i < ndofs; i++)
  {
    j = find_entry(Ai, i, Ap[i], Ap[i+1]);
    #ifndef COMPLEX
      if (j >= 0 && Ax[j] >= 0.0) equi[i] = sqrt(Ax[j]);
    #else
      if (j >= 0) equi[i] = sqrt(std::abs(Ax[j])); // todo: check
    #endif
  }

  // rescale the stiffness matrix and the RHS
  for (i = 0; i < ndofs; i++)
  {
    for (j = Ap[i]; j < Ap[i+1]; j++)
      Ax[j] /= equi[i] * equi[Ai[j]];
    RHS[i] /= equi[i];
  }

  is_equi = true;
}


void DiscreteProblem::precalc_equi_coefs()
{
  int i, m;
  memset(equi, 0, sizeof(double) * ndofs);
  verbose("Precalculating equilibration coefficients...");

  RefMap refmap;
  AsmList al;
  Element* e;

  for (m = 0; m < neq; m++)
  {
    PrecalcShapeset* fu = pss[m];
    BiForm* bf = biform[m] + m;
    Mesh* mesh = spaces[m]->get_mesh();

    for_all_active_elements(e, mesh)
    {
      update_limit_table(e->get_mode());
      fu->set_active_element(e);
      refmap.set_active_element(e);

      spaces[m]->get_element_assembly_list(e, &al);
      for (i = 0; i < al.cnt; i++)
      {
        if (al.dof[i] < 0) continue;
        fu->set_active_shape(al.idx[i]);
        scalar sy = 0.0, un = 0.0;
        if (bf->unsym) un = bf->unsym(fu, fu, &refmap, &refmap);
        if (bf->sym)   sy = bf->sym  (fu, fu, &refmap, &refmap);
        #ifndef COMPLEX
        equi[al.dof[i]] += (sy + un) * sqr(al.coef[i]);
        #else
        equi[al.dof[i]] += 0;//std::norm(sy + un) * sqr(al.coef[i]);
        #endif
      }
    }
  }

  for (i = 0; i < ndofs; i++)
    equi[i] = sqrt(equi[i]);
  is_equi = true;
}


//// solve /////////////////////////////////////////////////////////////////////////////////////////

bool DiscreteProblem::solve_system(int n, ...)
{
  double control[UMFPACK_CONTROL];
  umfpack_defaults(control);
  if (is_equi) control[UMFPACK_SCALE] = 0;

  if (!quiet) { verbose("Performing UMFPACK numeric analysis..."); begin_time(); }
  if (Numeric != NULL) umfpack_free_numeric(&Numeric);
  int status = ndofs ? umfpack_numeric(Ap, Ai, Ax, Symbolic, &Numeric, control, NULL) : UMFPACK_OK;
  if (!quiet) verbose("  (time: %g sec)", end_time());
  if (status != UMFPACK_OK)
  {
    umfpack_status(status);
    if (status != UMFPACK_WARNING_singular_matrix)
      return false;
  }

  free_solution_vector();
  vec = new scalar[ndofs];
  if (!quiet) { verbose("Solving the system..."); begin_time(); }
  status = ndofs ? umfpack_solve(UMFPACK_A, Ap, Ai, Ax, vec, RHS, Numeric, control, NULL) : UMFPACK_OK;
  if (!quiet) verbose("  (time: %g sec)", end_time());
  // todo: free Ax here?

  if (status != UMFPACK_OK)
  {
    umfpack_status(status);
    if (status != UMFPACK_WARNING_singular_matrix)
      return false;
  }

  if (is_equi)
    for (int i = 0; i < ndofs; i++)
      vec[i+1] /= equi[i];

  va_list ap;
  va_start(ap, n);
  if (n > neq) n = neq;
  for (int i = 0; i < n; i++)
  {
    Solution* sln = va_arg(ap, Solution*);
    sln->set_fe_solution(spaces[i], pss[i], vec);
  }
  va_end(ap);
  return true;
}


bool DiscreteProblem::solve_again(int n, ...)
{
  if (!quiet) { verbose("Solving the system..."); begin_time(); }
  int status = ndofs ? umfpack_solve(UMFPACK_A, Ap, Ai, Ax, vec+1, RHS, Numeric, NULL, NULL) : UMFPACK_OK;
  if (!quiet) verbose("  (time: %g sec)", end_time());
  vec[0] = 1.0; // "dirichlet dof"

  if (status != UMFPACK_OK)
  {
    umfpack_status(status);
    if (status != UMFPACK_WARNING_singular_matrix)
      return false;
  }

  if (is_equi)
    for (int i = 0; i < ndofs; i++)
      vec[i+1] /= equi[i];

  va_list ap;
  va_start(ap, n);
  if (n > neq) n = neq;
  for (int i = 0; i < n; i++)
  {
    Solution* sln = va_arg(ap, Solution*);
    sln->set_fe_solution(spaces[i], pss[i], vec);
  }
  va_end(ap);
  return true;
}


//// matrix output /////////////////////////////////////////////////////////////////////////////////

void DiscreteProblem::save_matrix_matlab(const char* filename, const char* varname)
{
  if (Ai == NULL || Ax == NULL) error("Matrix has not been created and/or assembled yet.");
  FILE* f = fopen(filename, "w");
  if (f == NULL) error("Could not open file %s for writing.", filename);
  verbose("Saving stiffness matrix in MATLAB format...");
  fprintf(f, "%% Size: %dx%d\n%% Nonzeros: %d\ntemp = zeros(%d, 3);\ntemp = [\n", ndofs, ndofs, Ap[ndofs], Ap[ndofs]);
  for (int j = 0; j < ndofs; j++)
    for (int i = Ap[j]; i < Ap[j+1]; i++)
      #ifndef COMPLEX
        fprintf(f, "%d %d %.18e\n", Ai[i]+1, j+1, Ax[i]);
      #else
        fprintf(f, "%d %d %.18e + %.18ei\n", Ai[i]+1, j+1, Ax[i].real(), Ax[i].imag());
      #endif
  fprintf(f, "];\n%s = spconvert(temp);\n", varname);
  fclose(f);
}

/* Saves the matrix in a coordinate format. */
void DiscreteProblem::save_matrix_coo(const char* filename)
{
  if (Ai == NULL || Ax == NULL) error("Matrix has not been created and/or assembled yet.");
  FILE* f = fopen(filename, "w");
  if (f == NULL) error("Could not open file %s for writing.", filename);
  verbose("Saving stiffness matrix in COO format...");
  fprintf(f, "%d %d %d\n", ndofs, ndofs, Ap[ndofs]);
  for (int j = 0; j < ndofs; j++)
    for (int i = Ap[j]; i < Ap[j+1]; i++)
      #ifndef COMPLEX
        fprintf(f, "%d %d %.18e\n", Ai[i]+1, j+1, Ax[i]);
      #else
        fprintf(f, "%d %d %.18e + %.18ei\n", Ai[i]+1, j+1, Ax[i].real(), Ax[i].imag());
      #endif
  fclose(f);
}


void DiscreteProblem::save_rhs_matlab(const char* filename, const char* varname)
{
  if (RHS == NULL) error("RHS has not been assembled yet.");
  FILE* f = fopen(filename, "w");
  if (f == NULL) error("Could not open file %s for writing.", filename);
  verbose("Saving RHS vector in MATLAB format...");
  fprintf(f, "%% Size: %dx1\n%s = [\n", ndofs, varname);
  for (int i = 0; i < ndofs; i++)
    #ifndef COMPLEX
      fprintf(f, "%.18e\n", RHS[i]);
    #else
      fprintf(f, "%.18e + %.18ei\n", RHS[i].real(), RHS[i].imag());
    #endif
  fprintf(f, "];\n");
  fclose(f);
}


void DiscreteProblem::save_matrix_bin(const char* filename)
{
  if (Ai == NULL || Ax == NULL) error("Matrix has not been created and/or assembled yet.");
  FILE* f = fopen(filename, "wb");
  if (f == NULL) error("Could not open file %s for writing.", filename);
  verbose("Saving stiffness matrix in binary format...");
  hermes2d_fwrite("H2DX\001\000\000\000", 1, 8, f);
  int ssize = sizeof(scalar), nnz = Ap[ndofs];
  hermes2d_fwrite(&ssize, sizeof(int), 1, f);
  hermes2d_fwrite(&ndofs, sizeof(int), 1, f);
  hermes2d_fwrite(&nnz, sizeof(int), 1, f);
  hermes2d_fwrite(Ap, sizeof(int), ndofs+1, f);
  hermes2d_fwrite(Ai, sizeof(int), nnz, f);
  hermes2d_fwrite(Ax, sizeof(scalar), nnz, f);
  fclose(f);
}


void DiscreteProblem::save_rhs_bin(const char* filename)
{
  if (RHS == NULL) error("RHS has not been assembled yet.");
  FILE* f = fopen(filename, "wb");
  if (f == NULL) error("Could not open file %s for writing.", filename);
  verbose("Saving RHS vector in binary format...");
  hermes2d_fwrite("H2DR\001\000\000\000", 1, 8, f);
  int ssize = sizeof(scalar);
  hermes2d_fwrite(&ssize, sizeof(int), 1, f);
  hermes2d_fwrite(&ndofs, sizeof(int), 1, f);
  hermes2d_fwrite(RHS, sizeof(scalar), ndofs, f);
  fclose(f);
}


//// interface /////////////////////////////////////////////////////////////////////////////////////

DiscreteProblem::DiscreteProblem()
{
  neq = 0;
  biform = NULL; liform = NULL;
  spaces = NULL; pss = NULL;
  Ai = Ap = NULL;
  Ax = NULL;
  RHS = NULL;
  Dir = NULL;
  Symbolic = NULL;
  Numeric = NULL;
  vec = NULL;
  extern_fns = NULL;
  num_extern = num_user_pss = 0;
  quiet = false;
  equi = NULL;
  is_equi = false;
}


void DiscreteProblem::set_num_equations(int neq)
{
  if (neq <= 0) error("Invalid number of equations.");
  if (neq > 10) warn("Large number of equations (%d). Is this the intent?", neq);
  free();

  // initialize the bilinear form
  biform = new BiForm*[neq];
  for (int i = 0; i < neq; i++)
  {
    biform[i] = new BiForm[neq];
    memset(biform[i], 0, sizeof(BiForm) * neq);
  }

  // init the rest of the arrays
  liform = new LiForm[neq];
  spaces = new Space*[neq];
  pss = new PrecalcShapeset*[neq];

  memset(liform, 0, sizeof(LiForm) * neq);
  memset(spaces, 0, sizeof(Space*) * neq);
  memset(pss, 0, sizeof(PrecalcShapeset*) * neq);

  this->neq = neq;
}


void DiscreteProblem::set_spaces(int n, ...)
{
  if (n <= 0 || n > neq) error("Bad number of spaces.");
  va_list ap;
  va_start(ap, n);
  for (int i = 0; i < neq; i++)
    spaces[i] = (i < n) ? va_arg(ap, Space*) : spaces[n-1];
  va_end(ap);
}


void DiscreteProblem::set_spaces(Space** spaces)
{
  memcpy(this->spaces, spaces, sizeof(Space*) * neq);
}


void DiscreteProblem::set_pss(int n, ...)
{
  if (n <= 0 || n > neq) error("Bad number of pss's.");

  va_list ap;
  va_start(ap, n);
  for (int i = 0; i < n; i++)
    pss[i] = va_arg(ap, PrecalcShapeset*);
  va_end(ap);
  num_user_pss = n;

  for (int i = n; i < neq; i++)
  {
    if (spaces[i]->get_shapeset() != spaces[n-1]->get_shapeset())
      error("Spaces with different shapesets must have different pss's.");
    pss[i] = new PrecalcShapeset(pss[n-1]);
  }
}


void DiscreteProblem::set_external_fns(int n, ...)
{
  num_extern = n;
  if (extern_fns != NULL) delete [] extern_fns;
  extern_fns = new MeshFunction*[n];

  va_list ap;
  va_start(ap, n);
  for (int i = 0; i < n; i++)
    extern_fns[i] = va_arg(ap, MeshFunction*);
  va_end(ap);
}


void DiscreteProblem::set_bilinear_form(int i, int j, BiFormFnVol unsym, BiFormFnVol sym, BiFormFnSurf surf)
{
  if (i < 0 || i >= neq || j < 0 || j >= neq)
    error("Bad equation number.");
  if ((unsym == BF_SYM || unsym == BF_ANTISYM) && (sym != NULL || surf != NULL))
    error("Illegal use of BF_SYM or BF_UNSYM");

  biform[i][j].unsym = unsym;
  biform[i][j].sym   = sym;
  biform[i][j].surf  = surf;
}


void DiscreteProblem::set_linear_form(int i, LiFormFnVol linear_form, LiFormFnSurf linear_form_surf)
{
  if (i < 0 || i >= neq) error("Bad equation number.");
  liform[i].lf   = linear_form;
  liform[i].surf = linear_form_surf;
}


void DiscreteProblem::copy(DiscreteProblem* ep)
{
  neq = ep->neq;

  // copy bilinear forms
  biform = new BiForm*[neq];
  for (int i = 0; i < neq; i++)
  {
    biform[i] = new BiForm[neq];
    memcpy(biform[i], ep->biform[i], sizeof(BiForm) * neq);
  }

  // copy the rest of the arrays
  liform = new LiForm[neq];
  spaces = new Space*[neq];
  pss = new PrecalcShapeset*[neq];

  memcpy(liform, ep->liform, sizeof(LiForm) * neq);
  memcpy(spaces, ep->spaces, sizeof(Space*) * neq);
  memcpy(pss, ep->pss, sizeof(PrecalcShapeset*) * neq);

  // copy external functions
  num_extern = ep->num_extern;
  extern_fns = new MeshFunction*[num_extern];
  memcpy(extern_fns, ep->extern_fns, sizeof(MeshFunction*) * num_extern);

  num_user_pss = neq;
  quiet = ep->quiet;
}


void DiscreteProblem::free_matrix_indices()
{
  if (Ap != NULL) { ::free(Ap); Ap = NULL; }
  if (Ai != NULL) { ::free(Ai); Ai = NULL; }
  if (Symbolic != NULL) { umfpack_free_symbolic(&Symbolic); Symbolic = NULL; }
  if (equi != NULL) { ::free(equi); equi = NULL; }
}

void DiscreteProblem::free_matrix_values()
{
  if (Ax != NULL)  { ::free(Ax);  Ax = NULL;  }
  if (RHS != NULL) { ::free(RHS); RHS = NULL; }
  if (Dir != NULL) { ::free(Dir); Dir = NULL; }
  if (Numeric != NULL) { umfpack_free_numeric(&Numeric); Numeric = NULL; }
}

void DiscreteProblem::free_solution_vector()
{
  if (vec != NULL) { delete [] vec; vec = NULL; }
}


void DiscreteProblem::free()
{
  free_matrix_indices();
  free_matrix_values();
  free_solution_vector();

  for (int i = 0; i < neq; i++)
    delete [] biform[i];

  delete [] biform;  biform = NULL;
  delete [] liform;  liform = NULL;
  delete [] spaces;  spaces = NULL;
  delete [] extern_fns;  extern_fns = NULL;

  for (int i = num_user_pss; i < neq; i++)
    delete pss[i];

  delete [] pss;  pss = NULL;

  neq = num_user_pss = num_extern = 0;
}


//// order limitation and warning //////////////////////////////////////////////////////////////////

/*void set_order_limit_table(int* tri_table, int* quad_table, int n)
{
  if (n < 24) error("Order limit tables must have at least 24 entries.");
  g_order_table_tri  = tri_table;
  g_order_table_quad = quad_table;
}


void update_limit_table(int mode)
{
  g_quad_2d_std.set_mode(mode);
  g_max_order = g_quad_2d_std.get_max_order();
  g_order_table = (mode == MODE_TRIANGLE) ? g_order_table_tri : g_order_table_quad;
}


void warn_order()
{
  if (!warned_order && warn_integration)
  {
    warn("Not enough integration rules for exact integration.");
    warned_order = true;
  }
}
*/
