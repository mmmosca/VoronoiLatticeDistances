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
