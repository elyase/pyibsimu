/*! \file timer.cpp
 *  \brief %Timer for cputime and realtime
 */

/* Copyright (c) 2005-2012 Taneli Kalvas. All rights reserved.
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

#include "timer.hpp"
#include "comptime.hpp"


Timer::Timer()
{
    _cpu1 = clock();
    ibs_gettimeofday( &_time1, NULL );
}


void Timer::start( void )
{
    _cpu1 = clock();
    ibs_gettimeofday( &_time1, NULL );
}


void Timer::stop( void )
{
    _cpu2 = clock();
    ibs_gettimeofday( &_time2, NULL );
}


double Timer::get_real_time( void ) const
{
    return( (double)(_time2.tv_sec - 
                     _time1.tv_sec) + 
            ((double)(_time2.tv_usec - 
                      _time1.tv_usec))/1.0e6 );
}


double Timer::get_cpu_time( void ) const
{
    return( (_cpu2-_cpu1)/(double)CLOCKS_PER_SEC );
}


std::ostream &operator<<( std::ostream &os, const Timer &t )
{
    os << t.get_cpu_time() << " s (" << t.get_real_time() << " s realtime)";
    return( os );
}