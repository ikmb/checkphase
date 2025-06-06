/*
 *    Copyright (C) 2018-2021 by Lars Wienbrandt and Jan Christian Kässens,
 *    Institute of Clinical Molecular Biology, Kiel University
 *    
 *    This file is part of EagleImp.
 *
 *    EagleImp is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    EagleImp is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with EagleImp. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <cstring>
#include <map>

extern "C" {
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <pwd.h>
#include <sys/inotify.h>
#include <unistd.h>
}

using namespace std;

template<typename T>
inline T roundToMultiple(T n, T mult) {
    if (n % mult)
        return n + mult - (n % mult);
    else
        return n;
}

template<typename T>
inline T reduceToMultiple(T n, T mult) {
    return n - (n % mult);
}

template<typename T>
inline T divideRounded(T a, T b) {
    return (a + b - 1) / b;
}

// convert NULL terminated C string to a reverse complemented std::string
// only charactes ATCG are supported, others are kept as they are
inline string reverseComplement(const char *s) {
    string rc(s);
    auto rcit = rc.begin();
    for (int i = rc.size()-1; i >= 0; i--) {
        char c = s[i];
        switch (c) {
        case 'A':
            c = 'T';
            break;
        case 'T':
            c = 'A';
            break;
        case 'C':
            c = 'G';
            break;
        case 'G':
            c = 'C';
            break;
        default:
            break;
        }
        *rcit = c;
        rcit++;
    }
    return rc;
}

#endif // UTILS_H
