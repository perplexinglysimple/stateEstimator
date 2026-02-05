/**
 * @file utils.h
 * @brief Lightweight logging and portability helpers.
 */
#ifndef UTILS_H
#define UTILS_H

#include <string.h>
#include <stdio.h>

// Single printf path avoids variadic dispatch issues across compilers.
#define LOG_MSG_PRINTF(level, ...) \
    do { \
        printf("[%s:%d] %s: ", __FILENAME__, __LINE__, level); \
        printf(__VA_ARGS__); \
        printf("\n"); \
    } while (0)

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

/** @brief Log the current function name at INFO level. */
#define LOG_FUNCTION() LOG_MSG_PRINTF("INFO", "%s called", __func__)

/** @brief Log an error message. */
#define LOG_ERROR(...) LOG_MSG_PRINTF("ERROR", __VA_ARGS__)
/** @brief Log a warning message. */
#define LOG_WARNING(...) LOG_MSG_PRINTF("WARNING", __VA_ARGS__)
/** @brief Log an info message. */
#define LOG_INFO(...) LOG_MSG_PRINTF("INFO", __VA_ARGS__)

#endif
