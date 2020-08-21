/*! \file polysolver.hpp
 *  \brief Polynomial solver
 */

/* Copyright (c) 2005-2010,2012 Taneli Kalvas. All rights reserved.
 *
 * You can redistribute this software and/or modify it under the terms
 * of the GNU General Public License as published by the Free Software
 * Foundation; either version 2 of the License, or (at your option)
 * any later version.
 * 
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this library (file "COPYING" included in the package);
 * if not, write to the Free Software Foundation, Inc., 51 Franklin
 * Street, Fifth Floor, Boston, MA 02110-1301 USA
 * 
 * If you have questions about your rights to use or distribute this
 * software, please contact Berkeley Lab's Technology Transfer
 * Department at TTD@lbl.gov. Other questions, comments and bug
 * reports should be sent directly to the author via email at
 * taneli.kalvas@jyu.fi.
 * 
 * NOTICE. This software was developed under partial funding from the
 * U.S.  Department of Energy.  As such, the U.S. Government has been
 * granted for itself and others acting on its behalf a paid-up,
 * nonexclusive, irrevocable, worldwide license in the Software to
 * reproduce, prepare derivative works, and perform publicly and
 * display publicly.  Beginning five (5) years after the date
 * permission to assert copyright is obtained from the U.S. Department
 * of Energy, and subject to any subsequent five (5) year renewals,
 * the U.S. Government is granted for itself and others acting on its
 * behalf a paid-up, nonexclusive, irrevocable, worldwide license in
 * the Software to reproduce, prepare derivative works, distribute
 * copies to the public, perform publicly and display publicly, and to
 * permit others to do so.
 */

#ifndef POLYSOLVER_HPP
#define POLYSOLVER_HPP 1


#include <stdint.h>


/*! \brief Solve quadric equation a*x^2 + b*x + c = 0.
 *
 *  Solves the quadric equation. Can also handle linear case if a is
 *  zero. Returns the number of roots found. Roots are filled in
 *  ascending order.
 */
uint32_t solve_quadratic( double a, double b, double c, 
			  double *x0, double *x1 );


/*! \brief Solve cubic equation a*x^3 + b*x^2 + c*x + d = 0.
 *
 *  Solves the cubic equation.  Can also handle lower order
 *  polynomials if coefficients are zero. Returns the number of
 *  roots found. Roots are filled in ascending order.
 */
uint32_t solve_cubic( double a, double b, double c, double d, 
		      double *x0, double *x1, double *x2 );


/*! \brief Solve quartic equation x^4 + a*x^3 + b*x^2 + c*x + d = 0.
 *
 *  Solves the quartic equation. Returns the number of roots
 *  found. Roots are filled in ascending order.
 */
uint32_t solve_quartic( double a, double b, double c, double d,
			double *x0, double *x1, double *x2, double *x3 );


#endif
