/**
 * @file utils.h
 * @brief Lightweight logging and portability helpers.
 */
#ifndef UTILS_H
#define UTILS_H

#include <string.h>
#include <stdio.h>

/*
 * Compile-time logging control for EKF/matrix library internals.
 * Default is enabled for easier bring-up/debug.
 * Set -DEKF_ENABLE_LOGS=1 in compiler flags to enable logging.
 */
#ifndef EKF_ENABLE_LOGS
#define EKF_ENABLE_LOGS 0
#endif

#if EKF_ENABLE_LOGS
#define LOG_MSG_PRINTF(level, ...)                                                                                     \
    do                                                                                                                 \
    {                                                                                                                  \
        printf("[%s:%d] %s: ", __FILENAME__, __LINE__, level);                                                         \
        printf(__VA_ARGS__);                                                                                           \
        printf("\n");                                                                                                  \
    } while (0)
#else
#define LOG_MSG_PRINTF(level, ...)                                                                                     \
    do                                                                                                                 \
    {                                                                                                                  \
    } while (0)
#endif

#ifndef __FILENAME__
#ifdef __linux__
#include <libgen.h>
#define GET_FILENAME(path) basename(path)
#elif defined(__TI_COMPILER_VERSION__)
#define GET_FILENAME(path) (path)
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
