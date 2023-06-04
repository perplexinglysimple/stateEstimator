#ifndef UTILS_H
#define UTILS_H

#include "string.h"
#include "stdio.h"

#define LOG_MSG_PRINTF(msg, level, ...) \
    do { \
        printf("[%s:%d] %s: ", __FILE__, __LINE__, level); \
        printf(msg __VA_OPT__(,) __VA_ARGS__); \
        printf("\n"); \
    } while (0)


#define LOG_ERROR(msg, ...) LOG_MSG_PRINTF(msg, "ERROR", __VA_ARGS__)
#define LOG_WARNING(msg, ...) LOG_MSG_PRINTF(msg, "WARNING", __VA_ARGS__)
#define LOG_INFO(msg, ...) LOG_MSG_PRINTF(msg, "INFO", __VA_ARGS__)

#endif