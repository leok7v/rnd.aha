#pragma once
#ifndef rt_header_included
#define rt_header_included

// nano runtime to make debugging, life, universe and everything a bit easier

#include <locale.h>
#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#ifdef _MSC_VER // cl.exe compiler:
#pragma warning(disable: 5045) // Compiler will insert Spectre mitigation
#pragma warning(disable: 4710) // '...': function not inlined
#pragma warning(disable: 4711) // function '...' selected for automatic inline expansion
#pragma warning(disable: 4820) // '...' bytes padding added after data member '...'
#endif

#define null ((void*)0) // like null_ptr better than NULL (0)

#define rt_countof(a) (sizeof(a) / sizeof((a)[0]))

#include "rt_generics.h"

typedef struct rt_debug_output_s {
    char  buffer[8 * 1024];
    char* rd; // read pointer
    char* wr; // write pointer
    int32_t max_prefix_len;
    int32_t max_function_len;
} rt_debug_output_t;

static void rt_output_line(const char* s) {
    const char* text = s;
    const char* c = strstr(text, "):");
    if (c != null) {
        while (c > text && *c != '\\') { c--; }
        if (c != text) { text = c + 1; }
    }
    static bool setlocale_called;
    if (!setlocale_called) { setlocale(LC_ALL, "en_US.UTF-8"); }
    #ifdef OutputDebugString // will be defined if Window.h header is included
        if (!setlocale_called) { SetConsoleOutputCP(CP_UTF8); }
        WCHAR utf16[4096];
        int n = MultiByteToWideChar(CP_UTF8, 0, s, -1,
                                    utf16, rt_countof(utf16) - 1);
        if (n > 0) {
            utf16[rt_countof(utf16) - 1] = 0x00;
            OutputDebugStringW(utf16);
        } else {
            OutputDebugStringA("UTF-8 to UTF-16 conversion error\n");
        }
        fprintf(stderr, "%s", text);
    #else
        fprintf(stderr, "%s", text);
    #endif
    setlocale_called = true;
}

static void rt_flush_buffer(rt_debug_output_t* out, const char* file,
                            int32_t line, const char* function) {
    if (out->wr > out->rd) {
        if ((out->wr - out->rd) >= (sizeof(out->buffer) - 4)) {
            strcpy(out->wr - 3, "...\n");
        }
        char prefix[1024];
        snprintf(prefix, sizeof(prefix) - 1, "%s(%d):", file, line);
        prefix[sizeof(prefix) - 1] = 0x00;
        char* start = out->rd;
        char* end = strchr(start, '\n');
        while (end != null) {
            *end = '\0';
            char output[8 * 1024];
            out->max_prefix_len = rt_max(out->max_prefix_len,
                                        (int32_t)strlen(prefix));
            out->max_function_len = rt_max(out->max_function_len,
                                           (int32_t)strlen(function));
            snprintf(output, sizeof(output) - 1, "%-*s %-*s %s\n",
                     (unsigned int)out->max_prefix_len, prefix,
                     (unsigned int)out->max_function_len, function,
                     start);
            output[sizeof(output) - 1] = 0x00;
            rt_output_line(output);
            start = end + 1;
            end = strchr(start, '\n');
        }
        // Move any leftover text to the beginning of the buffer
        size_t leftover_len = strlen(start);
        memmove(out->buffer, start, leftover_len + 1);
        out->rd = out->buffer;
        out->wr = out->buffer + leftover_len;
    }
}

static int32_t rt_vprintf_implementation(const char* file, int32_t line,
                                         const char* function, bool nl,
                                         const char* format, va_list args) {
    static thread_local rt_debug_output_t out;
    if (out.rd == null || out.wr == null) {
        out.rd = out.buffer;
        out.wr = out.buffer;
    }
    size_t space = sizeof(out.buffer) - (out.wr - out.buffer) - 4;
    int32_t n = vsnprintf(out.wr, space, format, args);
    if (n < 0) {
        rt_output_line("printf format error\n");
    } else {
        char* p = out.wr + n - 1;
        if (nl) {
            if (n == 0 || *p != '\n') { p++; *p++ = '\n'; *p++ = '\0'; n++; }
            rt_flush_buffer(&out, file, line, function);
        }
        if (n >= (int32_t)space) {
            // Handle buffer overflow
            strcpy(out.buffer + sizeof(out.buffer) - 4, "...");
            rt_flush_buffer(&out, file, line, function);
        } else {
            out.wr += n;
            if (strchr(out.wr - n, '\n')) {
                rt_flush_buffer(&out, file, line, function);
            }
        }
    }
    return n;
}

static int32_t rt_printf_implementation(const char* file, int32_t line,
                                        const char* func, bool nl,
                                        const char* format,
                                        ...) {
    va_list args;
    va_start(args, format);
    int32_t result = rt_vprintf_implementation(file, line, func, nl,
                                               format, args);
    va_end(args);
    return result;
}

#define rt_printf(...) rt_printf_implementation(__FILE__,           \
                       __LINE__, __func__, false, "" __VA_ARGS__)

#define rt_println(...) rt_printf_implementation(__FILE__,          \
                        __LINE__, __func__, true, "" __VA_ARGS__)

/*
    printf("\xF0\x9F\x98\x80 Hello\xF0\x9F\x91\x8B "
           "world\xF0\x9F\x8C\x8D!\n\xF0\x9F\x98\xA1 Goodbye "
           "\xF0\x9F\x98\x88 cruel \xF0\x9F\x98\xB1 "
           "Universe \xF0\x9F\x8C\xA0\xF0\x9F\x8C\x8C..."
           "\xF0\x9F\x92\xA4\n");

    rt_printf("\xF0\x9F\x98\x80 Hello\xF0\x9F\x91\x8B ");
    rt_printf("world\xF0\x9F\x8C\x8D!\n\xF0\x9F\x98\xA1 Goodbye ");
    rt_printf("\xF0\x9F\x98\x88 cruel \xF0\x9F\x98\xB1 ");
    rt_printf("Universe \xF0\x9F\x8C\xA0\xF0\x9F\x8C\x8C...");
    rt_printf("\xF0\x9F\x92\xA4\n");

    rt_printf("\xF0\x9F\x98\x80 Hello\xF0\x9F\x91\x8B "
              "world\xF0\x9F\x8C\x8D!\n\xF0\x9F\x98\xA1 Goodbye "
              "\xF0\x9F\x98\x88 cruel \xF0\x9F\x98\xB1 "
              "Universe \xF0\x9F\x8C\xA0\xF0\x9F\x8C\x8C..."
              "\xF0\x9F\x92\xA4\n");
*/

#if defined(_MSC_VER)
    #define rt_swear(b, ...) ((void)                                        \
    ((!!(b)) || rt_printf_implementation(__FILE__, __LINE__, __func__,      \
                         true, #b " false " __VA_ARGS__) &&                 \
                rt_exit(1)))
#else
    #define rt_swear(b, ...) ((void)                                        \
        ((!!(b)) || rt_printf_implementation(__FILE__, __LINE__, __func__,  \
                         true, #b " false " __VA_ARGS__) &&                 \
                    rt_exit(1)))
#endif

#if defined(DEBUG) || defined(_DEBUG)
#define rt_assert(b, ...) rt_swear(b, __VA_ARGS__)
#else
#define rt_assert(b, ...) ((void)(0))
#endif

static int32_t rt_exit(int exit_code) {
    _Pragma("warning(push)")
    _Pragma("warning(disable: 4702)") /* unreachable code */
    if (exit_code == 0) { rt_printf("exit code must not be zero"); }
    if (exit_code != 0) {
        #ifdef _WINDOWS_
            DebugBreak();
            ExitProcess(exit_code);
        #else
            exit(exit_code);
        #endif
    }
    return 0;
    _Pragma("warning(pop)")
}

#endif // rt_header_included
