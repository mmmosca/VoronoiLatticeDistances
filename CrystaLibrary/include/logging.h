/*Copyright (C) 2026 Marco M. Mosca

This file is part of VoronoiLatticeDistances.

VoronoiLatticeDistances is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

VoronoiLatticeDistances is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with VoronoiLatticeDistances. If not, see <https://www.gnu.org/licenses/>.
*/
#ifndef _LOGGING_H
#define _LOGGING_H

#include <string>

enum class LogLevel {
    LOGDEBUG = 0,
    LOGINFO = 1,
    LOGWARNING = 2,
    LOGERROR = 3,
    LOGFATAL = 4
};

class Logger
{
public:
    static void init(const std::string& filename);

    static void setLevel(LogLevel inputLevel);
    static void debug(const std::string& message);
    static void info(const std::string& message);
    static void warn(const std::string& message);
    static void error(const std::string& message);

private:
    static void log(const std::string& level,
                    const std::string& message);
};

#endif // !_LOGGING_H
