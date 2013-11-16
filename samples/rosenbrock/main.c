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

// This program shows how to use nmead.h to minimize one version of the 
// Rosenbrock function with 3 decision variables.

#include <stdio.h>
#include "nmead.h"

//Rosenbrock function with 3 variables. It has a single minimum at (1, 1, 1) for
//any c > 0.
//Notice that it accepts a pointer to a double (array with values of independent
//variables) and a void pointer. You can use the latter to pass extra
//parameters to the function.
double f(const double * restrict x, void * restrict par)
{
  //Get the constant.
  const double c = *((double *)par);

  return (1 - x[0]) * (1 - x[0]) + 
            c * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]) +
              (1 - x[1]) * (1 - x[1]) +
                  c * (x[2] - x[1] * x[1]) * (x[2] - x[1] * x[1]);

}

//"Instantiate" necessary functions and objects for libnmead. We use double and
//3 decision variables.  
NMEAD_DEFINE_PROTOTYPES(double, 3)

int main()
{
  //Initial guess.
  double x0[] = {0.5, -0.5, 0.5};

  //Tolerance level. This controls when to stop.
  const double tol = 1e-7;

  //"Step size". This controls how we build the initial simplex.
  const double step = 0.1;

  //Rosenbrock function parameter (the traditional value is 1.0, but the
  //solution will be the same for any c > 0).
  double c = 3.0;

  //Get a minimizer using a macro. Notice that this is very similar to the use
  //of a template object on C++, where the "template parameters" are the type of
  //the involved scalars, and the number of decision variables.
  NMEAD_SOLVER(double, 3) s;

  //Initialize the minimizer.
  NMEAD_SOLVER_INIT(double, 3)(&s, &f, &c, x0, step);
    
  //Iterate. 
  size_t it = 0;
  //Maximum number of iterations.
  const double maxit = 500;

  while (it <= maxit) {
    //Update the guess.
    NMEAD_SOLVER_STEP(double, 3)(&s);

    //Print current guess, function value, and simplex "size" (the smaller, the
    //better).
    const double *xmin = NMEAD_SOLVER_MIN(double, 3)(&s);
    const double fmin = NMEAD_SOLVER_FMIN(double, 3)(&s);
    
    printf("it = %lu, x = [%f, %f, %f], f = %f, size = %f\n", 
            it, xmin[0], xmin[1], xmin[2], fmin, 
            NMEAD_SOLVER_SIZE(double, 3)(&s)); 
  
    //Check for convergence.
    if (NMEAD_SOLVER_DONE(double, 3)(&s, tol)) { break; }

    ++it;
  } 
}
