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
#include "logging.h"

#include <fstream>
#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <type_traits>

namespace
{
    std::unique_ptr<std::ofstream> logFile;
    std::ostream* logHandle = &std::cerr;
    LogLevel level = LogLevel::LOGERROR;
    unsigned short int levelValue = static_cast<unsigned short int>(level);
}

void Logger::init(const std::string& filename)
{
    if (!logFile->is_open())
    {
        logFile = std::make_unique<std::ofstream>(
            filename,
            std::ios::app);
        logHandle = logFile.get();
    }
}

void Logger::setLevel(LogLevel inputLevel) {
    level = inputLevel;
    levelValue = static_cast<unsigned short int>(level);
}


void Logger::debug(const std::string& message)
{
    if (static_cast<unsigned short int>(LogLevel::LOGDEBUG) >= levelValue) {
        log("DEBUG", message);
    }
}

void Logger::info(const std::string& message)
{
    if (static_cast<unsigned short int>(LogLevel::LOGINFO) >= levelValue) {
        log("INFO", message);
    }
}

void Logger::warn(const std::string& message)
{
    if (static_cast<unsigned short int>(LogLevel::LOGWARNING) >= levelValue) {
        log("WARNING", message);
    }
}

void Logger::error(const std::string& message)
{
    if (static_cast<unsigned short int>(LogLevel::LOGERROR) >= levelValue) {
        log("ERROR", message);
    }
}

void Logger::log(const std::string& level,
                 const std::string& message)
{
    auto now = std::chrono::system_clock::now();
    auto time = std::chrono::system_clock::to_time_t(now);

    (*logHandle) << std::put_time(std::localtime(&time),
                               "%Y-%m-%d %H:%M:%S")
              << " [" << level << "] "
              << message
              << '\n';
}
