#pragma once
    #ifdef _DEBUG_
    #   define ASSERT_MSG(condition, message) \
        do { \
            if (! (condition)) { \
                std::cerr << "Assertion `" #condition "` failed in " << __FILE__ \
                          << " line " << __LINE__ << ": \nError!!! " << message << std::endl; \
                std::terminate(); \
            } \
        } while (false)
    #else
    #   define ASSERT_MSG(condition, message) do { } while (false)
    #endif


    #ifdef _DEBUG_
    #define DEBUG_MSG(str) do { std::cout << str << std::endl; } while( false )
    #else
    #define DEBUG_MSG(str) do { } while ( false )
    #endif    
