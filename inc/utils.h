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

#ifdef __FILENAME__ 
#define LOG_MSG_PRINTF(msg, level, ...) \
    do { \
        printf("[%s:%d] %s: ", __FILENAME__, __LINE__, level); \
        _VA_OPT_PRINTF_CHECK(msg, __VA_ARGS__); \
        printf("\n"); \
    } while (0)
#else
#define LOG_MSG_PRINTF(msg, level, ...) \
    do { \
        printf("[%s:%d] %s: ", __FILE__, __LINE__, level); \
        _VA_OPT_PRINTF_CHECK(msg, __VA_ARGS__); \
        printf("\n"); \
    } while (0)
#endif


#define LOG_ERROR(msg, ...) LOG_MSG_PRINTF(msg, "ERROR", __VA_ARGS__)
#define LOG_WARNING(msg, ...) LOG_MSG_PRINTF(msg, "WARNING", __VA_ARGS__)
#define LOG_INFO(msg, ...) LOG_MSG_PRINTF(msg, "INFO", __VA_ARGS__)

#endif