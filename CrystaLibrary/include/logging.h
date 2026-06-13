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
