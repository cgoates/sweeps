#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <set>
#include <vector>
#include <chrono>
#include <csignal>

// clang-format off
#define LOG( COND ) if( COND ) std::cout
// clang-format on

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

class Timer
{
    public:
    Timer( const size_t n = 10 )
        : mAccumulated( n, std::chrono::microseconds::zero() ),
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
            mAccumulated.at( i ) += std::chrono::duration_cast<std::chrono::microseconds>(
                std::chrono::high_resolution_clock::now() - mStartTimes.at( i ) );
            mStarted.at( i ) = false;
        }
        return mAccumulated.at( i ).count() * 1e-6;
    }

    private:
    std::vector<std::chrono::microseconds> mAccumulated;
    std::vector<std::chrono::time_point<std::chrono::high_resolution_clock>> mStartTimes;
    std::vector<bool> mStarted;
};

bool equals( const double& a, const double& b, const double& tol );

bool equals( const Eigen::Ref<const Eigen::VectorXd> a, const Eigen::Ref<const Eigen::VectorXd> b, const double& tol );