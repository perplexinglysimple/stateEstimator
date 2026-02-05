/**
 * @file utils.h
 * @brief Lightweight logging and portability helpers.
 */
#ifndef UTILS_H
#define UTILS_H

#include "string.h"
#include "stdio.h"

#define VA_OPT_SUPPORTED

#if __cplusplus <= 201703 && defined __GNUC__ \
  && !defined __clang__ && !defined __EDG__ // These compilers pretend to be GCC
#undef VA_OPT_SUPPORTED
#endif

#ifndef VL_STUDIO_STUPIDITY
#undef VA_OPT_SUPPORTED
#endif


#if defined(VA_OPT_SUPPORTED)
#define _VA_OPT_PRINTF_CHECK(msg, ...) \
    printf(msg __VA_OPT__(,) __VA_ARGS__)
#else
#define _VA_OPT_PRINTF_CHECK(msg, ...) \
    printf(msg, ##__VA_ARGS__)
#endif

#define PRINTF_1ARG(msg) \
    printf(msg)

#define PRINTF_2ARGS(msg, ...) \
    printf(msg, __VA_ARGS__)

#define GET_3RD_ARG(arg1, arg2, arg3, ...) arg3

#define _LOG_MSG_MACRO_CHOOSER(...) \
    GET_3RD_ARG(__VA_ARGS__, PRINTF_1ARGS, PRINTF_2ARG, )( __VA_ARGS__)

#ifndef __FILENAME__
#ifdef __linux__
#include <libgen.h>
#define GET_FILENAME(path) basename(path)
#elif _WIN32
#include <windows.h>
#include <tchar.h>
#define MAX_PATH_LENGTH 256

/** @brief Extract the filename portion from a Windows path. */
const char* GetFileNameFromPath(const char* path);

#define GET_FILENAME(path) GetFileNameFromPath(path)
#else
#error "Unsupported platform"
#endif
#define __FILENAME__ GET_FILENAME(__FILE__)
#endif

#define LOG_MSG_PRINTF(msg, level, ...) \
    do { \
        printf("[%s:%d] %s: ", __FILENAME__, __LINE__, level); \
        _VA_OPT_PRINTF_CHECK(msg, __VA_ARGS__); \
        printf("\n"); \
    } while (0)

/** @brief Log the current function name at INFO level. */
#define LOG_FUNCTION() LOG_MSG_PRINTF("%s called", "INFO", __func__)

/** @brief Log an error message. */
#define LOG_ERROR(msg, ...) LOG_MSG_PRINTF(msg, "ERROR", __VA_ARGS__)
/** @brief Log a warning message. */
#define LOG_WARNING(msg, ...) LOG_MSG_PRINTF(msg, "WARNING", __VA_ARGS__)
/** @brief Log an info message. */
#define LOG_INFO(msg, ...) LOG_MSG_PRINTF(msg, "INFO", __VA_ARGS__)

#endif
