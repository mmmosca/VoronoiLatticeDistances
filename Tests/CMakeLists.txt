include(FetchContent)

FetchContent_Declare(
    googletest
    URL https://github.com/google/googletest/releases/download/v1.17.0/googletest-1.17.0.tar.gz
)

# Prevent overriding parent compiler settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(googletest)

add_executable(
    tests
    test_crystalstructure.cpp
    test_voronoidiagram.cpp
)

target_link_libraries(
    tests
    PRIVATE
    CrystaLibrary
    GTest::gtest_main
)

gtest_discover_tests(tests)
