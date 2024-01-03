#pragma once
#include<iostream>
#include<Eigen/Dense>
#include<Eigen/Sparse>
#include<set>
#include<vector>

#define LOG( COND ) if( COND ) std::cout


std::ostream& operator<<( std::ostream& o, const Eigen::Triplet<double>& t );

template<typename T>
std::ostream& operator<<( std::ostream& o, const std::vector<T>& v )
{
    if( v.size() == 0 ) o << "{}";
    else
    {    
        o << "{ ";
        for( auto it = v.begin(); it != v.end() - 1; it++ ) o << *it << ", ";
        o << v.back() << " }";
    }
    return o;
}

std::ostream& operator<<( std::ostream& o, const std::vector<Eigen::Vector3d>& v );

template<typename T>
std::ostream& operator<<( std::ostream& o, const std::set<T>& v )
{
    if( v.size() == 0 ) o << "{}";
    else
    {
        o << "{ ";
        for( auto it = v.begin(); it != v.end(); it++ ) o << *it << ", ";
        o << " }";
    }
    return o;
}