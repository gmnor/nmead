// nmead.h

// Copyright (C) 2013 Gilberto Noronha
//
// This is free software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.

// This software is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.

// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 51
// Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

#ifndef NMEAD_H_
#define NMEAD_H_

//For bool.
#include <stdbool.h>

//For size_t.
#include <stdio.h>

//Type name.
#define NMEAD_POINT(T, N)  nmead_point_ ## T ## _ ## N 

//Minimizer.
#define NMEAD_SOLVER(T, N)  nmead_solver_ ## T ## _ ## N 

//Function names for internal stuff. 
#define NMEAD_F(T, N, FUNC) nmead_point_ ## T ## _ ## N ## _ ## FUNC

//Function names for interface stuff.
#define NMEAD_FUNC(T, N, FUNC) nmead_ ## FUNC ## _ ## T ## _ ## N     

//Comparator.
#define NMEAD_CMP(T, N) nmead_point_ ## T ## _ ## N ## _cmp

//We need to implement abs because libc uses different names for different
//types.
#define NMEAD_ABS(T) nmead_abs_ ## T

//Convenience macros to make it easier to call the more important functions.
#define NMEAD_COMPUTE_SIZE(T, N) NMEAD_FUNC(T, N, compute_size)
#define NMEAD_STEP(T, N) NMEAD_FUNC(T, N, step)
#define NMEAD_SOLVE(T, N) NMEAD_FUNC(T, N, solve)

#define NMEAD_SORT(T, N) NMEAD_F(T, N, qsort)

#define NMEAD_SOLVER_INIT(T, N) NMEAD_FUNC(T, N, solver_init)
#define NMEAD_SOLVER_STEP(T, N) NMEAD_FUNC(T, N, solver_step)
#define NMEAD_SOLVER_MIN(T, N) NMEAD_FUNC(T, N, solver_min)
#define NMEAD_SOLVER_FMIN(T, N) NMEAD_FUNC(T, N, solver_fmin)
#define NMEAD_SOLVER_DONE(T, N) NMEAD_FUNC(T, N, solver_done)
#define NMEAD_SOLVER_SIZE(T, N) NMEAD_FUNC(T, N, solver_size)

#define NMEAD_DEFINE_PROTOTYPES(T, N) \
  \
  /*Generic implementation of abs.*/\
  static inline T NMEAD_ABS(T)(const T x)\
  {\
    return x < (T)0 ? -x : x;\
  }\
  \
  /*Helper type*/\
  typedef struct\
  {\
    T data[N + 1];\
  } NMEAD_POINT(T, N);\
  \
  \
  /*Copy points.*/\
  static inline void NMEAD_F(T, N, copy)(const NMEAD_POINT(T, N) * restrict src,\
                                         NMEAD_POINT(T, N) * restrict dest)\
  {\
    for (size_t i = 0; i < N + 1; ++i) {\
      dest->data[i] = src->data[i];\
    }\
  }\
  \
  /* 
   * Weighted average of two points (z = x + c * (x - y)).
   * Notice that we DO NOT COPY THE LAST COORDINATE since that is the function
   * evaluation.
   * */\
  static inline void NMEAD_F(T, N, reflect)(\
                              const T c,\
                              const NMEAD_POINT(T, N) * restrict x,\
                              const NMEAD_POINT(T, N) * restrict y,\
                              NMEAD_POINT(T, N) * restrict z)\
  {\
    for (size_t i = 0; i < N; ++i) {\
      z->data[i] = x->data[i] + c * (x->data[i] - y->data[i]);\
    }\
  }\
  \
  /*Comparator for sorting.*/\
  static inline int NMEAD_F(T, N, cmp)(const void * x, \
                                       const void * y)\
  {\
    return ((NMEAD_POINT(T, N) *)(x))->data[N] >\
           ((NMEAD_POINT(T, N) *)(y))->data[N];\
  }\
  \
  /*Swap two points.*/\
  static inline void NMEAD_F(T, N, swap)(NMEAD_POINT(T, N) *a,\
                                    NMEAD_POINT(T, N) *b)\
  {\
    for (size_t i = 0; i < N + 1; ++i) {\
      const T tmp = a->data[i];\
      a->data[i] = b->data[i];\
      b->data[i] = tmp;\
    }\
  }\
  \
  /* Given two points, swap them if the first is greater than the second.*
   * TODO: It might be better to make this a macro.*/\
  static inline void NMEAD_F(T, N, compare_and_swap)(NMEAD_POINT(T, N) *a,\
                                    NMEAD_POINT(T, N) *b)\
  {\
    if (a->data[N] > b->data[N]) {\
      NMEAD_F(T, N, swap)(a, b);\
    }\
  }\
  \
  /* Sort 3 points using network sort.
   * TODO: Check whether other combinations are faster.
  static inline void NMEAD_F(T, 2, sort_network)(NMEAD_POINT(T, 2) *m)\
  {\
    NMEAD_F(T, N, compare_and_swap)(m, m + 1);\
    NMEAD_F(T, N, compare_and_swap)(m, m + 2);\
    NMEAD_F(T, N, compare_and_swap)(m + 1, m + 2);\
  }\*/\
  \
  /* Sort 4 points using network sort.
   * TODO: Check whether other combinations are faster.*/\
  static inline void NMEAD_F(T, 3, sort_network)(NMEAD_POINT(T, 3) *m)\
  {\
    NMEAD_F(T, N, compare_and_swap)(m, m + 1);\
    NMEAD_F(T, N, compare_and_swap)(m + 2, m + 3);\
    NMEAD_F(T, N, compare_and_swap)(m, m + 2);\
    NMEAD_F(T, N, compare_and_swap)(m + 1, m + 3);\
    NMEAD_F(T, N, compare_and_swap)(m + 1, m + 2);\
  }\
  \
  static inline size_t NMEAD_F(T, N, qsort_helper)(NMEAD_POINT(T, N) *first,\
                                                   const size_t right,\
                                                   const size_t ind)\
  {\
    /*Copy current pivot.*/\
    NMEAD_POINT(T, N) tmp;\
    for (size_t i = 0; i < N + 1; ++i) {\
      tmp.data[i] = (first + ind)->data[i];\
    }\
    \
    /*Move pivot to the end of the partition.*/\
    NMEAD_F(T, N, swap)(first + ind, first + right);\
    \
    /*Move elements smaller than the pivot to the left.*/\
    size_t store = 0;\
    \
    for (size_t i = 0; i < right; ++i) {\
      /*if ((first + i)->data[N] < tmp.data[N]) {\*/\
      if (!NMEAD_F(T, N, cmp)(first + i, &tmp)) {\
        NMEAD_F(T, N, swap)(first + i, first + store);\
        ++store;\
      }\
    }\
    \
    /*Move pivot to its correct position.*/\
    NMEAD_F(T, N, swap)(first + store, first + right);\
    \
    /*Return pivot's position.*/\
    return store;\
  }\
  \
  static inline void NMEAD_F(T, N, qsort)(NMEAD_POINT(T, N) *m, size_t n)\
  {\
    /*We need at least three elements to use recursive calls.*/\
    if (n >= 5) {\
      /*
       *  Partition the sequence. Notice that we use use the middle point (n / 2 -
       *  1) as the pivot. This can be optimized...
       *  
       *  */\
      size_t ind = NMEAD_F(T, N, qsort_helper)(m, n - 1,  n / 2 - 1);\
      \
      /*Sort the two partitions.*/\
      NMEAD_F(T, N, qsort)(m, ind + 1);\
      NMEAD_F(T, N, qsort)(m + ind + 1, n - ind - 1);\
    }\
    \
    else if (n == 4) {\
      NMEAD_F(T, 3, sort_network)(m);\
    }\
    /*If there are only two elements left we compare them.*/\
    else if (n == 2) {\
      NMEAD_F(T, N, compare_and_swap)(m, m + 1);\
    }\
    /*With just one element, there is nothing to do.*/\
    else {\
      return;\
    }\
  }\
  \
  \
  /*Iteration function.*/\
  \
  static void NMEAD_STEP(T, N)(T (*f)(const T *, void *),\
                               NMEAD_POINT(T, N) *data,\
                               void *params)\
  {\
    /*Coefficients of reflection, expansion, contraction and shrink,
     * respectively.*/\
    const T alpha = 1.0;\
    const T gamma = 2.0;\
    const T rho = -0.5;\
    const T sigma = 0.5;\
    \
    \
    /*Sort the points.*/\
    /*NMEAD_F(T, N, sort_network)(data);\*/\
    NMEAD_F(T, N, qsort)(data, N + 1);\
    /*qsort(data, N + 1, sizeof(NMEAD_POINT(T, N)), &NMEAD_F(T, N, cmp));*/\
    \
    /*Compute center of simplex.*/\
    NMEAD_POINT(T, N) x0;\
    for (size_t i = 0; i < N; ++i) {\
      x0.data[i] = data[0].data[i];\
      for (size_t j = 1; j < N; ++j) {\
        x0.data[i] += data[j].data[i];\
      }\
      x0.data[i] /= ((T)N);\
    }\
    /*
     * Reflection. We want to see if a point far from the center is useful.
     * 
    */\
    NMEAD_POINT(T, N) xr;\
    /*xr = (1 + alpha) * x0 - alpha * data[N].*/\
    NMEAD_F(T, N, reflect)(alpha, &x0, &data[N], &xr);\
    \
    /*Evaluate the function at the reflected point.*/\
    xr.data[N] = f(xr.data, params);\
    /*Replace the worst point with the reflected if the latter is better than
     * the second worst, but not better than the best.*/\
    if (xr.data[N] < data[N - 1].data[N] && xr.data[N] >= data[0].data[N]) {\
      NMEAD_F(T, N, copy)(&xr, &data[N]);\
    }\
    \
    /*If the reflected point is better than the best, we expand the simplex.*/\
    else if (xr.data[N] < data[0].data[N]) {\
      NMEAD_POINT(T, N) xe;\
      /*xe = (1 + gamma) * x0 - gamma * data[N].*/\
      NMEAD_F(T, N, reflect)(gamma, &x0, &data[N], &xe);\
      \
      /*Evaluate the function a second time.*/\
      xe.data[N] = f(xe.data, params);\
      \
      /*If the expanded is better than the reflected, we replace the worst with
       * the expanded. Otherwise, we replace the worst with the reflected.*/\
      if (xe.data[N] < xr.data[N]) {\
        NMEAD_F(T, N, copy)(&xe, &data[N]);\
      }\
      else {\
        NMEAD_F(T, N, copy)(&xr, &data[N]);\
      }\
    }\
    \
    /*If we get here, we know that the reflected point is not better than the
     * second-worst. In this case, we try contraction.*/\
    else {\
      NMEAD_POINT(T, N) xc;\
      /*xc = (1 + rho) * x0 - rho * data[N].*/\
      NMEAD_F(T, N, reflect)(rho, &x0, &data[N], &xc);\
      \
      /*Evaluate the function a second time.*/\
      xc.data[N] = f(xc.data, params);\
      \
      /*If the contracted point is better than the worst, we replace the worst 
        * with the contracted.*/\
      if (xc.data[N] < data[N].data[N]) {\
        NMEAD_F(T, N, copy)(&xc, &data[N]);\
      }\
      \
      /*Bad (and uncommon) case: The contracted point is not the better than any
       * current point. In this case, we replace all points except for the best.
       */\
      else {\
        for (size_t i = 0; i < N; ++i) {\
          NMEAD_F(T, N, reflect)(-sigma, &data[0], &data[i + 1], &data[i + 1]);\
          /*Need to evaluate the function N times!*/\
          data[i + 1].data[N] = f(data[i + 1].data, params);\
        }\
      }\
    }\
    \
  }\
  /*
   * Helper function to compute the "size" of the simplex.
  */\
  static inline T NMEAD_COMPUTE_SIZE(T, N)(NMEAD_POINT(T, N) *points)\
  {\
    T size = 0;\
    for (size_t i = 0; i < N; ++i) {\
      const T tmp = points[0].data[i];\
      for (size_t j = 1; j < N + 1; ++j) {\
        size += NMEAD_ABS(T)(points[j].data[i] - tmp);\
      }\
    }\
    \
    return size / ((T)(N * N));\
  }\
  \
  /*
   * Main minimization function.
  */\
  static inline void NMEAD_SOLVE(T, N)(T *guesses, void *params,\
                                     T (*f)(const T *, void*),\
                                     const size_t maxit, const T tol)\
  {\
    /*Copy guesses to array of points.*/\
    NMEAD_POINT(T, N) data[N + 1];\
    for (size_t i = 0; i < N + 1; ++i) {\
      for (size_t j = 0; j < N; ++j) {\
        /*Assume that the guesses are stored in row-major order.*/\
        data[i].data[j] = guesses[i * N + j];\
      }\
      /*Need to evaluate the function. at each guess point.*/\
      data[i].data[N] = f(data[i].data, params);\
    }\
    \
    /*Now we iterate.*/\
    size_t it = 0;\
    while (it <= maxit) {\
      NMEAD_STEP(T, N)(f, data, params);\
      \
      /*Check if we can stop.*/\
      if (NMEAD_COMPUTE_SIZE(T, N)(data) < tol) { break; }\
      \
      ++it;\
   }\
   for (size_t i = 0; i < N; ++i) {\
    guesses[i] = data[0].data[i];\
   }\
 }\
 \
 /*Minimizer object.*/\
 typedef struct\
 {\
   /*Simplex points.*/\
   NMEAD_POINT(T, N) x[N + 1];\
   \
   /*"Size" of simplex.*/\
   T size;\
   \
   /*Function to be minimized.*/\
   T (*f)(const T *, void *);\
   \
   /*Parameters to be passed to the function.*/\
   void *params;\
  } NMEAD_SOLVER(T, N);\
  \
  \
  /*The next functions provide the public API to use solvers to minimize
   * functions.*/\
  static inline void NMEAD_FUNC(T, N, solver_init)(\
                                          NMEAD_SOLVER(T, N) *self,\
                                          T (*f)(const T *, void *),\
                                          void *params,\
                                          const T * guess,\
                                          const T step_size)\
  {\
    /*Set function and parameters.*/\
    self->f = f;\
    self->params = params;\
    \
    /*Copy guesses.*/\
    for (size_t i = 0; i < N + 1; ++i) {\
      for (size_t j = 0; j < N; ++j) {\
        self->x[i].data[j] = guess[j];\
      }\
    }\
    /*Evaluate function at the first point.*/\
    self->x[0].data[N] = f(self->x[0].data, params);\
    \
    /*The other points are adjusted using the provided step_size.*/\
    for (size_t i = 0; i < N; ++i) {\
      self->x[i + 1].data[i] += step_size;\
      /*Need to evaluate the function. at each guess point.*/\
      self->x[i + 1].data[N] = f(self->x[i + 1].data, params);\
    }\
    /*Initialize the size to a big number, but do NOT compute it (only do it
     * when the user explicitly asks us to do so)*/\
    self->size = (T)(1000000);\
  }\
  \
  /*Compute current simplex size.*/\
  static inline T NMEAD_FUNC(T, N, solver_size)(\
                                          NMEAD_SOLVER(T, N) *self)\
  {\
    self->size = NMEAD_COMPUTE_SIZE(T, N)(self->x);\
    return self->size;\
  }\
  /*Return current minimum estimate.*/\
  static inline T* NMEAD_FUNC(T, N, solver_min)(\
                            NMEAD_SOLVER(T, N) *self)\
  {\
    return self->x[0].data;\
  }\
  /*Return function value at the minimum.*/\
  static inline T NMEAD_FUNC(T, N, solver_fmin)(\
                            NMEAD_SOLVER(T, N) *self)\
  {\
    return self->x[0].data[N];\
  }\
  /*Iterate once.*/\
  static inline void NMEAD_FUNC(T, N, solver_step)(\
                                          NMEAD_SOLVER(T, N) *self)\
  {\
    NMEAD_STEP(T, N)(self->f, self->x, self->params);\
  }\
  /*Test for convergence.*/\
  static inline bool NMEAD_FUNC(T, N, solver_done)(\
                                          NMEAD_SOLVER(T, N) *self,\
                                          const T tol)\
  {\
    return (tol >= NMEAD_FUNC(T, N, solver_size)(self));\
  }\
  \
  /*Iterate at most maxit times.*/\
  static inline int NMEAD_FUNC(T, N, solver_stepmany)(\
                                          NMEAD_SOLVER(T, N) *self,\
                                          const T tol,\
                                          const size_t maxit)\
  {\
    size_t it = 0;\
    while (!NMEAD_FUNC(T, N, solver_done)(self, tol) && it <= maxit) {\
      NMEAD_FUNC(T, N, solver_step)(self);\
      ++it;\
    }\
    return (it < maxit);\
  }

#endif // NMEAD_H_
