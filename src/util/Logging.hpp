#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <set>
#include <vector>
#include <chrono>
#include <optional>
#include <SmallVector.hpp>

// clang-format off
#define LOG( COND ) if( COND ) std::cout
// clang-format on

namespace util
{
    inline std::string black() { return "\033[30m"; }
    inline std::string red() { return "\033[31m"; }
    inline std::string green() { return "\033[32m"; }
    inline std::string yellow() { return "\033[33m"; }
    inline std::string blue() { return "\033[34m"; }
    inline std::string purple() { return "\033[35m"; }
    inline std::string cyan() { return "\033[36m"; }
    inline std::string bblack() { return "\033[1;30m"; }
    inline std::string bred() { return "\033[1;31m"; }
    inline std::string bgreen() { return "\033[1;32m"; }
    inline std::string byellow() { return "\033[1;33m"; }
    inline std::string bblue() { return "\033[1;34m"; }
    inline std::string bpurple() { return "\033[1;35m"; }
    inline std::string bcyan() { return "\033[1;36m"; }
    inline std::string reset() { return "\033[0m"; }
    inline std::string resetl() { return "\033[0m\n"; }
}

std::ostream& operator<<( std::ostream& o, const Eigen::Triplet<double>& t );

template <typename T> std::ostream& operator<<( std::ostream& o, const std::vector<T>& v )
{
    if( v.size() == 0 )
        o << "{}";
    else
    {
        o << "{ ";
        for( auto it = v.begin(); it != v.end() - 1; it++ ) o << *it << ", ";
        o << v.back() << " }";
    }
    return o;
}

std::ostream& operator<<( std::ostream& o, const std::vector<Eigen::Vector3d>& v );

template <typename T> std::ostream& operator<<( std::ostream& o, const std::set<T>& v )
{
    if( v.size() == 0 )
        o << "{}";
    else
    {
        o << "{ ";
        for( auto it = v.begin(); it != v.end(); it++ ) o << *it << ", ";
        o << " }";
    }
    return o;
}

template <typename T> std::ostream& operator<<( std::ostream& o, const std::optional<T>& v )
{
    if( not v.has_value() )
        o << "{}";
    else
    {
        o << v.value();
    }
    return o;
}

template <typename T, typename U> std::ostream& operator<<( std::ostream& o, const std::map<T, U>& v )
{
    if( v.size() == 0 )
        o << "{}";
    else
    {
        o << "{ ";
        for( auto it = v.begin(); it != v.end(); it++ ) o << "{ " << it->first << ", " << it->second << " }" << ", ";
        o << " }";
    }
    return o;
}

void pauseDebugger();

class Timer
{
    public:
    Timer( const size_t n = 10 )
        : mAccumulated( n, std::chrono::nanoseconds::zero() ),
          mStartTimes( n, std::chrono::high_resolution_clock::now() ),
          mStarted( n, false )
    {}

    void start( const size_t i )
    {
        mStartTimes.at( i ) = std::chrono::high_resolution_clock::now();
        mStarted.at( i ) = true;
    }

    double stop( const size_t i )
    {
        if( mStarted.at( i ) )
        {
            mAccumulated.at( i ) += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::high_resolution_clock::now() - mStartTimes.at( i ) );
            mStarted.at( i ) = false;
        }
        return mAccumulated.at( i ).count() * 1e-9;
    }

    private:
    std::vector<std::chrono::nanoseconds> mAccumulated;
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> mStartTimes;
    std::vector<bool> mStarted;
};