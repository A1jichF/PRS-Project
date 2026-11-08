#pragma once

// Minimalni test framework bez vanjskih zavisnosti.
//   TEST(ime) { CHECK(uslov); CHECK_NEAR(a, b, tol); }

#include <cmath>
#include <cstdio>
#include <functional>
#include <string>
#include <vector>

namespace testing {

struct TestCase {
    const char*           name;
    std::function<void()> fn;
};

inline std::vector<TestCase>& registry()
{
    static std::vector<TestCase> tests;
    return tests;
}

inline int& failureCount()
{
    static int count = 0;
    return count;
}

struct Registrar {
    Registrar(const char* name, std::function<void()> fn) { registry().push_back({name, std::move(fn)}); }
};

inline void reportFailure(const char* file, int line, const std::string& msg)
{
    ++failureCount();
    std::fprintf(stderr, "  FAIL %s:%d: %s\n", file, line, msg.c_str());
}

} // namespace testing

#define TEST_CONCAT_INNER(a, b) a##b
#define TEST_CONCAT(a, b) TEST_CONCAT_INNER(a, b)

#define TEST(name)                                                                      \
    static void name();                                                                 \
    static ::testing::Registrar TEST_CONCAT(registrar_, name)(#name, name);             \
    static void name()

#define CHECK(cond)                                                                     \
    do {                                                                                \
        if (!(cond))                                                                    \
            ::testing::reportFailure(__FILE__, __LINE__, "CHECK(" #cond ")");           \
    } while (0)

#define CHECK_NEAR(a, b, tol)                                                           \
    do {                                                                                \
        const double va_ = static_cast<double>(a);                                      \
        const double vb_ = static_cast<double>(b);                                      \
        if (!(std::abs(va_ - vb_) <= static_cast<double>(tol)))                         \
            ::testing::reportFailure(__FILE__, __LINE__,                                \
                std::string("CHECK_NEAR(" #a ", " #b ") : ") + std::to_string(va_) +     \
                " vs " + std::to_string(vb_));                                          \
    } while (0)
