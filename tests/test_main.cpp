#include "TestFramework.h"

#include <cstdio>
#include <exception>

int main()
{
    int failedTests = 0;
    for (const auto& t : testing::registry()) {
        const int before = testing::failureCount();
        try {
            t.fn();
        } catch (const std::exception& e) {
            testing::reportFailure(__FILE__, __LINE__, std::string("izuzetak: ") + e.what());
        }
        const bool ok = testing::failureCount() == before;
        if (!ok)
            ++failedTests;
        std::printf("[%s] %s\n", ok ? " OK " : "FAIL", t.name);
    }
    std::printf("\n%zu testova, %d palo\n", testing::registry().size(), failedTests);
    return failedTests == 0 ? 0 : 1;
}
