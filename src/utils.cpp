/* utils.cpp
 * 
 * Copyright (C) 2018  Marco van Zwetselaar <io@zwets.it>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <string>
#include <cstdio>
#include <cstdarg>
#include <cstdlib>
#include <thread>
#include <unistd.h>
#include "utils.h"

namespace kfc {

static bool verbose = false;
static const char* progname = "";

void
set_progname(const char *p)
{
    if (p)
        progname = p;
}

bool
set_verbose(bool v)
{
    bool old_verbose = verbose;
    verbose = v;
    return old_verbose;
}

void
emit(const char *fmt, ...)
{
    char buf[2048];

    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);

    std::cerr << progname << ": " << buf << std::endl;
}

void
raise_error(const char *fmt, ...)
{
    char buf[2048];

    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);

    std::cerr << progname << ": error: " << buf << std::endl;
    std::exit(1);
}

void
verbose_emit(const char *fmt, ...)
{
    if (verbose)
    {
        char buf[2048];

        va_list ap;
        va_start(ap, fmt);
        vsnprintf(buf, sizeof(buf), fmt, ap);
        va_end(ap);

        std::cerr << progname << ": " << buf << std::endl;
    }
}

// NOTE: get_system_memory may fail to compile on MacOS.
// TODO: solve that issue with ifdefs
//
unsigned long long
get_system_memory()
{
    return (static_cast<unsigned long long>(sysconf(_SC_PHYS_PAGES)) * static_cast<unsigned long long>(sysconf(_SC_PAGE_SIZE)));
}

unsigned int
get_system_threads()
{
    unsigned int nthreads = std::thread::hardware_concurrency();
    if (nthreads == 0) {
        std::cerr << progname << ": warning: cannot determine hardware concurrency, using 1 thread" << std::endl;
        nthreads = 1;
    }
    return nthreads;
}


} // namespace kfc

// vim: sts=4:sw=4:ai:si:et
