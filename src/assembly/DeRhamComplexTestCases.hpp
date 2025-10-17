#pragma once
#include <DeRhamComplexDiscretization.hpp>

using namespace assembly;

enum class TestCase
{
    Case2d,
    MoreThanMinimalIntersection_2d,
    MinimalIntersection_2d,
    LessThanMinimalIntersection_2d,
    ShortestChainExists_2d,
    Case2d_nForm,
    Case4x4xn_3d,
    Case2p1KittyCorners_3d,
    Case2p1KittyCornersIntersection_3d,
    Case2p1KittyCornersIntersection_3dfixed,
    Case2p1KittyCornersIntersection_2d
};

DeRhamComplexHierarchicalDiscretization cases( const TestCase which_case, std::optional<size_t> parameter = std::nullopt )
{
    switch( which_case )
    {
        case TestCase::Case2d:
        {
            const size_t degree = 4;

    // const basis::KnotVector kv_s( util::concatenate({ std::vector<double>( degree, 0.0 ),
    //                                                              util::linspace( 0.0, std::numbers::pi, 11 ),
    //                                                              std::vector<double>( degree, std::numbers::pi ) }), 1e-9 );
    const basis::KnotVector kv_s = basis::unitIntervalKnotVectorWithNElems( 16, degree );
    // const basis::KnotVector kv_s = basis::integerKnotsWithNElems( 2, degree );

    const basis::KnotVector kv_t = kv_s;
    const param::ParentDomain domain = param::cubeDomain( 2 );
    const assembly::QuadratureRule quad_rule = assembly::QuadratureRule( degree + 1, domain );
    const Eigen::MatrixXd control_points = util::tensorProduct( { grevillePoints( kv_s, degree ), grevillePoints( kv_t, degree ) } ).transpose();
    // DeRhamComplexTPDiscretization drcd( kv_s, kv_t, degree, degree, control_points );
    DeRhamComplexHierarchicalDiscretization drcd( kv_s,
                                                  kv_t,
                                                  degree,
                                                  degree,
                                                  control_points,
                                                //   { { { 1, 1 },
                                                //       { 1, 2 },
                                                //       { 1, 3 },
                                                //       { 1, 4 },
                                                //       { 1, 5 },
                                                //       { 1, 6 },
                                                //       { 1, 7 },
                                                //       { 2, 1 },
                                                //       { 2, 2 },
                                                //       { 2, 3 },
                                                //       { 2, 4 },
                                                //       { 2, 5 },
                                                //       { 2, 6 },
                                                //       { 2, 7 },
                                                //       { 3, 1 },
                                                //       { 3, 2 },
                                                //       { 3, 3 },
                                                //       { 3, 4 },
                                                //       { 3, 5 },
                                                //       { 3, 6 },
                                                //       { 3, 7 },
                                                //       { 4, 1 },
                                                //       { 4, 2 },
                                                //       { 4, 3 },
                                                //       { 4, 5 },
                                                //       { 4, 6 },
                                                //       { 4, 7 },
                                                //       { 5, 3 },
                                                //       { 5, 4 },
                                                //       { 5, 5 },
                                                //       { 5, 6 },
                                                //       { 5, 7 },
                                                //       { 6, 3 },
                                                //       { 6, 4 },
                                                //       { 6, 5 },
                                                //       { 6, 6 },
                                                //       { 6, 7 },
                                                //       { 7, 3 },
                                                //       { 7, 4 },
                                                //       { 7, 5 },
                                                //       { 7, 6 },
                                                //       { 7, 7 }
                                                //      } } );
    //   { { { 1, 1 },
    //       { 1, 2 },
    //       { 1, 3 },
    //       { 2, 1 },
    //       { 2, 2 },
    //       { 2, 3 },
    //       { 3, 1 },
    //       { 3, 2 },
    //       { 3, 3 },
    //       { 3, 4 },
    //       { 3, 5 },
    //       { 4, 3 },
    //       { 4, 4 },
    //       { 4, 5 },
    //       { 5, 3 },
    //       { 5, 4 },
    //       { 5, 5 } } } );
    { {
        { 5, 0 }, { 5, 1 }, { 5, 2 }, { 5, 3 }, { 5, 4 }, { 5, 5 }, { 5, 6 }, { 5, 7 }, { 5, 8 }, { 5, 9 }, { 5, 10 }, { 5, 11 }, { 5, 12 }, { 5, 13 }, { 5, 14 }, { 5, 15 },
        { 6, 0 }, { 6, 1 }, { 6, 2 }, { 6, 3 }, { 6, 4 }, { 6, 5 }, { 6, 6 }, { 6, 7 }, { 6, 8 }, { 6, 9 }, { 6, 10 }, { 6, 11 }, { 6, 12 }, { 6, 13 }, { 6, 14 }, { 6, 15 },
        { 7, 0 }, { 7, 1 }, { 7, 2 }, { 7, 3 }, { 7, 4 }, { 7, 5 }, { 7, 6 }, { 7, 7 }, { 7, 8 }, { 7, 9 }, { 7, 10 }, { 7, 11 }, { 7, 12 }, { 7, 13 }, { 7, 14 }, { 7, 15 }
    }

    } );
            return drcd;
        }
        case TestCase::MoreThanMinimalIntersection_2d:
        {
            const size_t degree_s = 3;
            const size_t degree_t = 4;
            const size_t n_elems_1d_s = 11;
            const size_t n_elems_1d_t = 11;
            const basis::KnotVector kv_s = basis::unitIntervalKnotVectorWithNElems( n_elems_1d_s, degree_s );
            const basis::KnotVector kv_t = basis::unitIntervalKnotVectorWithNElems( n_elems_1d_t, degree_t );
            const param::ParentDomain domain = param::cubeDomain( 2 );
            const assembly::QuadratureRule quad_rule = assembly::QuadratureRule( std::max( degree_s, degree_t ) + 1, domain );
            const Eigen::MatrixXd control_points = util::tensorProduct( { grevillePoints( kv_s, degree_s ), grevillePoints( kv_t, degree_t ) } ).transpose();

            std::vector<std::vector<std::pair<size_t, size_t>>> elems_to_refine;
            elems_to_refine.emplace_back();
            elems_to_refine.back().push_back( { 1, 1 } );
            elems_to_refine.back().push_back( { 1, 2 } );
            elems_to_refine.back().push_back( { 1, 3 } );
            elems_to_refine.back().push_back( { 1, 4 } );
            elems_to_refine.back().push_back( { 1, 5 } );

            elems_to_refine.back().push_back( { 2, 1 } );
            elems_to_refine.back().push_back( { 2, 2 } );
            elems_to_refine.back().push_back( { 2, 3 } );
            elems_to_refine.back().push_back( { 2, 4 } );
            elems_to_refine.back().push_back( { 2, 5 } );

            elems_to_refine.back().push_back( { 3, 1 } );
            elems_to_refine.back().push_back( { 3, 2 } );
            elems_to_refine.back().push_back( { 3, 3 } );
            elems_to_refine.back().push_back( { 3, 4 } );
            elems_to_refine.back().push_back( { 3, 5 } );

            elems_to_refine.back().push_back( { 4, 1 } );
            elems_to_refine.back().push_back( { 4, 2 } );
            elems_to_refine.back().push_back( { 4, 3 } );
            elems_to_refine.back().push_back( { 4, 4 } );
            elems_to_refine.back().push_back( { 4, 5 } );

            elems_to_refine.back().push_back( { 4, 8 } );
            elems_to_refine.back().push_back( { 4, 7 } );
            elems_to_refine.back().push_back( { 4, 6 } );
            elems_to_refine.back().push_back( { 3, 6 } );
            elems_to_refine.back().push_back( { 3, 7 } );
            elems_to_refine.back().push_back( { 3, 8 } );

            elems_to_refine.back().push_back( { 5, 8 } );
            elems_to_refine.back().push_back( { 5, 4 } );
            elems_to_refine.back().push_back( { 5, 5 } );
            elems_to_refine.back().push_back( { 5, 6 } );
            elems_to_refine.back().push_back( { 5, 7 } );

            elems_to_refine.back().push_back( { 6, 8 } );
            elems_to_refine.back().push_back( { 6, 4 } );
            elems_to_refine.back().push_back( { 6, 5 } );
            elems_to_refine.back().push_back( { 6, 6 } );
            elems_to_refine.back().push_back( { 6, 7 } );

            return DeRhamComplexHierarchicalDiscretization( kv_s,
                                                        kv_t,
                                                        degree_s,
                                                        degree_t,
                                                        control_points,
                                                        elems_to_refine );
        }
        case TestCase::MinimalIntersection_2d:
        {
            const size_t degree_s = 3;
            const size_t degree_t = 4;
            const size_t n_elems_1d_s = 11;
            const size_t n_elems_1d_t = 11;
            const basis::KnotVector kv_s = basis::unitIntervalKnotVectorWithNElems( n_elems_1d_s, degree_s );
            const basis::KnotVector kv_t = basis::unitIntervalKnotVectorWithNElems( n_elems_1d_t, degree_t );
            const param::ParentDomain domain = param::cubeDomain( 2 );
            const assembly::QuadratureRule quad_rule = assembly::QuadratureRule( std::max( degree_s, degree_t ) + 1, domain );
            const Eigen::MatrixXd control_points = util::tensorProduct( { grevillePoints( kv_s, degree_s ), grevillePoints( kv_t, degree_t ) } ).transpose();

            std::vector<std::vector<std::pair<size_t, size_t>>> elems_to_refine;
            elems_to_refine.emplace_back();
            elems_to_refine.back().push_back( { 1, 1 } );
            elems_to_refine.back().push_back( { 1, 2 } );
            elems_to_refine.back().push_back( { 1, 3 } );
            elems_to_refine.back().push_back( { 1, 4 } );
            elems_to_refine.back().push_back( { 1, 5 } );

            elems_to_refine.back().push_back( { 2, 1 } );
            elems_to_refine.back().push_back( { 2, 2 } );
            elems_to_refine.back().push_back( { 2, 3 } );
            elems_to_refine.back().push_back( { 2, 4 } );
            elems_to_refine.back().push_back( { 2, 5 } );

            elems_to_refine.back().push_back( { 3, 1 } );
            elems_to_refine.back().push_back( { 3, 2 } );
            elems_to_refine.back().push_back( { 3, 3 } );
            elems_to_refine.back().push_back( { 3, 4 } );
            elems_to_refine.back().push_back( { 3, 5 } );

            elems_to_refine.back().push_back( { 4, 1 } );
            elems_to_refine.back().push_back( { 4, 2 } );
            elems_to_refine.back().push_back( { 4, 3 } );
            elems_to_refine.back().push_back( { 4, 4 } );
            elems_to_refine.back().push_back( { 4, 5 } );

            elems_to_refine.back().push_back( { 5, 8 } );
            elems_to_refine.back().push_back( { 5, 4 } );
            elems_to_refine.back().push_back( { 5, 5 } );
            elems_to_refine.back().push_back( { 5, 6 } );
            elems_to_refine.back().push_back( { 5, 7 } );

            elems_to_refine.back().push_back( { 6, 8 } );
            elems_to_refine.back().push_back( { 6, 4 } );
            elems_to_refine.back().push_back( { 6, 5 } );
            elems_to_refine.back().push_back( { 6, 6 } );
            elems_to_refine.back().push_back( { 6, 7 } );

            elems_to_refine.back().push_back( { 7, 8 } );
            elems_to_refine.back().push_back( { 7, 4 } );
            elems_to_refine.back().push_back( { 7, 5 } );
            elems_to_refine.back().push_back( { 7, 6 } );
            elems_to_refine.back().push_back( { 7, 7 } );

            elems_to_refine.back().push_back( { 8, 8 } );
            elems_to_refine.back().push_back( { 8, 4 } );
            elems_to_refine.back().push_back( { 8, 5 } );
            elems_to_refine.back().push_back( { 8, 6 } );
            elems_to_refine.back().push_back( { 8, 7 } );

            return DeRhamComplexHierarchicalDiscretization( kv_s,
                                                        kv_t,
                                                        degree_s,
                                                        degree_t,
                                                        control_points,
                                                        elems_to_refine );
        }
        case TestCase::LessThanMinimalIntersection_2d:
        {
            const size_t degree_s = 3;
            const size_t degree_t = 4;
            const size_t n_elems_1d_s = 11;
            const size_t n_elems_1d_t = 11;
            const basis::KnotVector kv_s = basis::unitIntervalKnotVectorWithNElems( n_elems_1d_s, degree_s );
            const basis::KnotVector kv_t = basis::unitIntervalKnotVectorWithNElems( n_elems_1d_t, degree_t );
            const param::ParentDomain domain = param::cubeDomain( 2 );
            const assembly::QuadratureRule quad_rule = assembly::QuadratureRule( std::max( degree_s, degree_t ) + 1, domain );
            const Eigen::MatrixXd control_points = util::tensorProduct( { grevillePoints( kv_s, degree_s ), grevillePoints( kv_t, degree_t ) } ).transpose();

            std::vector<std::vector<std::pair<size_t, size_t>>> elems_to_refine;
            elems_to_refine.emplace_back();
            elems_to_refine.back().push_back( { 1, 1 } );
            elems_to_refine.back().push_back( { 1, 2 } );
            elems_to_refine.back().push_back( { 1, 3 } );
            elems_to_refine.back().push_back( { 1, 4 } );
            elems_to_refine.back().push_back( { 1, 5 } );

            elems_to_refine.back().push_back( { 2, 1 } );
            elems_to_refine.back().push_back( { 2, 2 } );
            elems_to_refine.back().push_back( { 2, 3 } );
            elems_to_refine.back().push_back( { 2, 4 } );
            elems_to_refine.back().push_back( { 2, 5 } );

            elems_to_refine.back().push_back( { 3, 1 } );
            elems_to_refine.back().push_back( { 3, 2 } );
            elems_to_refine.back().push_back( { 3, 3 } );
            elems_to_refine.back().push_back( { 3, 4 } );
            elems_to_refine.back().push_back( { 3, 5 } );

            elems_to_refine.back().push_back( { 4, 1 } );
            elems_to_refine.back().push_back( { 4, 2 } );
            elems_to_refine.back().push_back( { 4, 3 } );
            elems_to_refine.back().push_back( { 4, 4 } );
            elems_to_refine.back().push_back( { 4, 5 } );

            elems_to_refine.back().push_back( { 5, 9 } );
            elems_to_refine.back().push_back( { 5, 5 } );
            elems_to_refine.back().push_back( { 5, 6 } );
            elems_to_refine.back().push_back( { 5, 7 } );
            elems_to_refine.back().push_back( { 5, 8 } );

            elems_to_refine.back().push_back( { 6, 9 } );
            elems_to_refine.back().push_back( { 6, 5 } );
            elems_to_refine.back().push_back( { 6, 6 } );
            elems_to_refine.back().push_back( { 6, 7 } );
            elems_to_refine.back().push_back( { 6, 8 } );

            elems_to_refine.back().push_back( { 7, 9 } );
            elems_to_refine.back().push_back( { 7, 5 } );
            elems_to_refine.back().push_back( { 7, 6 } );
            elems_to_refine.back().push_back( { 7, 7 } );
            elems_to_refine.back().push_back( { 7, 8 } );

            elems_to_refine.back().push_back( { 8, 9 } );
            elems_to_refine.back().push_back( { 8, 5 } );
            elems_to_refine.back().push_back( { 8, 6 } );
            elems_to_refine.back().push_back( { 8, 7 } );
            elems_to_refine.back().push_back( { 8, 8 } );

            return DeRhamComplexHierarchicalDiscretization( kv_s,
                                                        kv_t,
                                                        degree_s,
                                                        degree_t,
                                                        control_points,
                                                        elems_to_refine );
        }
        case TestCase::ShortestChainExists_2d:
        {
            const size_t degree_s = 3;
            const size_t degree_t = 4;
            const size_t n_elems_1d_s = 11;
            const size_t n_elems_1d_t = 11;
            const basis::KnotVector kv_s = basis::unitIntervalKnotVectorWithNElems( n_elems_1d_s, degree_s );
            const basis::KnotVector kv_t = basis::unitIntervalKnotVectorWithNElems( n_elems_1d_t, degree_t );
            const Eigen::MatrixXd control_points = util::tensorProduct( { grevillePoints( kv_s, degree_s ), grevillePoints( kv_t, degree_t ) } ).transpose();

            std::vector<std::vector<std::pair<size_t, size_t>>> elems_to_refine;
            elems_to_refine.emplace_back();
            elems_to_refine.back().push_back( { 1, 1 } );
            elems_to_refine.back().push_back( { 1, 2 } );
            elems_to_refine.back().push_back( { 1, 3 } );
            elems_to_refine.back().push_back( { 1, 4 } );
            elems_to_refine.back().push_back( { 1, 5 } );

            elems_to_refine.back().push_back( { 2, 1 } );
            elems_to_refine.back().push_back( { 2, 2 } );
            elems_to_refine.back().push_back( { 2, 3 } );
            elems_to_refine.back().push_back( { 2, 4 } );
            elems_to_refine.back().push_back( { 2, 5 } );

            elems_to_refine.back().push_back( { 3, 1 } );
            elems_to_refine.back().push_back( { 3, 2 } );
            elems_to_refine.back().push_back( { 3, 3 } );
            elems_to_refine.back().push_back( { 3, 4 } );
            elems_to_refine.back().push_back( { 3, 5 } );

            elems_to_refine.back().push_back( { 4, 1 } );
            elems_to_refine.back().push_back( { 4, 2 } );
            elems_to_refine.back().push_back( { 4, 3 } );
            elems_to_refine.back().push_back( { 4, 4 } );
            elems_to_refine.back().push_back( { 4, 5 } );

            elems_to_refine.back().push_back( { 4, 8 } );
            elems_to_refine.back().push_back( { 4, 7 } );
            elems_to_refine.back().push_back( { 4, 6 } );
            elems_to_refine.back().push_back( { 3, 6 } );
            elems_to_refine.back().push_back( { 3, 7 } );
            elems_to_refine.back().push_back( { 3, 8 } );

            elems_to_refine.back().push_back( { 5, 8 } );
            elems_to_refine.back().push_back( { 5, 4 } );
            elems_to_refine.back().push_back( { 5, 5 } );
            elems_to_refine.back().push_back( { 5, 6 } );
            elems_to_refine.back().push_back( { 5, 7 } );

            elems_to_refine.back().push_back( { 6, 8 } );
            elems_to_refine.back().push_back( { 6, 4 } );
            elems_to_refine.back().push_back( { 6, 5 } );
            elems_to_refine.back().push_back( { 6, 6 } );
            elems_to_refine.back().push_back( { 6, 7 } );

            elems_to_refine.back().push_back( { 5, 2 } );
            elems_to_refine.back().push_back( { 5, 3 } );
            elems_to_refine.back().push_back( { 6, 3 } );
            elems_to_refine.back().push_back( { 1, 6 } );
            elems_to_refine.back().push_back( { 2, 6 } );
            elems_to_refine.back().push_back( { 2, 7 } );


            return DeRhamComplexHierarchicalDiscretization( kv_s,
                                                        kv_t,
                                                        degree_s,
                                                        degree_t,
                                                        control_points,
                                                        elems_to_refine );
        }
        case TestCase::Case2d_nForm:
        {
            const size_t degree = 4;
            const size_t n_elems_1d = 8;
            const basis::KnotVector kv_s = basis::unitIntervalKnotVectorWithNElems( n_elems_1d, degree );

            const basis::KnotVector kv_t = kv_s;
            const param::ParentDomain domain = param::cubeDomain( 2 );
            const assembly::QuadratureRule quad_rule = assembly::QuadratureRule( degree + 1, domain );
            const Eigen::MatrixXd control_points = util::tensorProduct( { grevillePoints( kv_s, degree ), grevillePoints( kv_t, degree ) } ).transpose();

            std::vector<std::vector<std::pair<size_t, size_t>>> elems_to_refine;
            elems_to_refine.emplace_back();
            const size_t start_elem = ( n_elems_1d - degree ) / 2;
            const size_t end_elem = ( n_elems_1d + degree ) / 2;
            for( size_t j = start_elem; j < end_elem; j++ )
                for( size_t k = start_elem; k < end_elem; k++ )
                    elems_to_refine.back().push_back( { j, k } );

            if( parameter.has_value() )
            {
                for( size_t level = 1; level < parameter.value(); level++ )
                {
                    elems_to_refine.emplace_back();
                    const size_t level_num_elems_1d = n_elems_1d * pow( 2, level );
                    const size_t start_elem = ( level_num_elems_1d - degree ) / 2;
                    const size_t end_elem = ( level_num_elems_1d + degree ) / 2;
                    for( size_t j = start_elem; j < end_elem; j++ )
                        for( size_t k = start_elem; k < end_elem; k++ )
                            elems_to_refine.back().push_back( { j, k } );
                }
            }
            DeRhamComplexHierarchicalDiscretization drcd( kv_s,
                                                        kv_t,
                                                        degree,
                                                        degree,
                                                        control_points,
                                                        elems_to_refine );
            return drcd;
        }
        case TestCase::Case4x4xn_3d:
        {
            const size_t degree = 4;
            const size_t n_elems_1d = 8;
            const basis::KnotVector kv_s = basis::unitIntervalKnotVectorWithNElems( n_elems_1d, degree );

            const basis::KnotVector kv_t = kv_s;
            const basis::KnotVector kv_u = kv_s;
            const Eigen::MatrixXd control_points = util::tensorProduct( { grevillePoints( kv_s, degree ), grevillePoints( kv_t, degree ), grevillePoints( kv_u, degree ) } ).transpose();

            std::vector<std::vector<std::array<size_t, 3>>> elems_to_refine;
            elems_to_refine.emplace_back();
            for( size_t i = 0; i < n_elems_1d; i++ )
            {
                for( size_t j = n_elems_1d / 2 - 2; j <= n_elems_1d / 2 + 1; j++ )
                    for( size_t k = n_elems_1d / 2 - 2; k <= n_elems_1d / 2 + 1; k++ )
                        elems_to_refine.back().push_back( { j, k, i } );
            }

            if( parameter.has_value() )
            {
                for( size_t level = 2; level <= parameter.value(); level++ )
                {
                    elems_to_refine.emplace_back();
                    size_t level_num_elems_1d = n_elems_1d * level;
                    for( size_t i = 0; i < level_num_elems_1d; i++ )
                    {
                        for( size_t j = level_num_elems_1d / 2 - 2; j <= level_num_elems_1d / 2 + 1; j++ )
                            for( size_t k = level_num_elems_1d / 2 - 2; k <= level_num_elems_1d / 2 + 1; k++ )
                                elems_to_refine.back().push_back( { j, k, i } );
                    }
                }
            }

            DeRhamComplexHierarchicalDiscretization drcd( kv_s, kv_t, kv_u, degree, degree, degree, control_points, { elems_to_refine } );
            return drcd;
        }
        case TestCase::Case2p1KittyCorners_3d:
        {
            const size_t degree = 3;
            const size_t n_elems_1d = 2 * degree + 2;
            const basis::KnotVector kv_s = basis::unitIntervalKnotVectorWithNElems( n_elems_1d, degree );

            const basis::KnotVector kv_t = kv_s;
            const basis::KnotVector kv_u = kv_s;
            const Eigen::MatrixXd control_points = util::tensorProduct( { grevillePoints( kv_s, degree ), grevillePoints( kv_t, degree ), grevillePoints( kv_u, degree ) } ).transpose();

            std::vector<std::array<size_t, 3>> elems_to_refine;
            for( size_t i = 0; i <= degree; i++ )
            {
                for( size_t j = 0; j <= degree; j++ )
                {
                    for( size_t k = 0; k <= degree; k++ )
                    {
                        elems_to_refine.push_back( { j, k, i } );
                        elems_to_refine.push_back( { n_elems_1d - 1 - j, n_elems_1d - 1 - k, n_elems_1d - 1 - i } );
                    }
                }
            }
            std::cout << "Refining " << elems_to_refine << std::endl;

            DeRhamComplexHierarchicalDiscretization drcd( kv_s, kv_t, kv_u, degree, degree, degree, control_points, { elems_to_refine } );
            return drcd;
        }
        case TestCase::Case2p1KittyCornersIntersection_3d:
        {
            const size_t degree = 3;
            const size_t n_elems_1d = 2 * degree + 2;
            const basis::KnotVector kv_s = basis::unitIntervalKnotVectorWithNElems( n_elems_1d, degree );

            const basis::KnotVector kv_t = kv_s;
            const basis::KnotVector kv_u = kv_s;
            const Eigen::MatrixXd control_points = util::tensorProduct( { grevillePoints( kv_s, degree ), grevillePoints( kv_t, degree ), grevillePoints( kv_u, degree ) } ).transpose();

            std::vector<std::array<size_t, 3>> elems_to_refine;
            for( size_t i = 0; i <= degree; i++ )
            {
                for( size_t j = 0; j < 2*degree; j++ )
                {
                    for( size_t k = 0; k < 2*degree; k++ )
                    {
                        elems_to_refine.push_back( { j, k, i } );
                        if( j <= degree and k <= degree )
                            elems_to_refine.push_back( { n_elems_1d - 1 - j, n_elems_1d - 1 - k, n_elems_1d - 1 - i } );
                    }
                }
            }

            DeRhamComplexHierarchicalDiscretization drcd( kv_s, kv_t, kv_u, degree, degree, degree, control_points, { elems_to_refine } );
            return drcd;
        }
        case TestCase::Case2p1KittyCornersIntersection_3dfixed:
        {
            const size_t degree = 3;
            const size_t n_elems_1d = 2 * degree + 2;
            const basis::KnotVector kv_s = basis::unitIntervalKnotVectorWithNElems( n_elems_1d, degree );

            const basis::KnotVector kv_t = kv_s;
            const basis::KnotVector kv_u = kv_s;
            const Eigen::MatrixXd control_points = util::tensorProduct( { grevillePoints( kv_s, degree ), grevillePoints( kv_t, degree ), grevillePoints( kv_u, degree ) } ).transpose();

            std::vector<std::array<size_t, 3>> elems_to_refine;
            for( size_t i = 0; i <= degree; i++ )
            {
                for( size_t j = 0; j <= 2*degree; j++ )
                {
                    for( size_t k = 0; k <= 2*degree; k++ )
                    {
                        elems_to_refine.push_back( { j, k, i } );
                        if( ( j <= degree and k <= degree + 1 ) or ( k > 0 and k <= degree + 1 and j <= degree + 1 ) )
                            elems_to_refine.push_back( { n_elems_1d - 1 - j, n_elems_1d - 1 - k, n_elems_1d - 1 - i } );
                    }
                }
            }

            DeRhamComplexHierarchicalDiscretization drcd( kv_s, kv_t, kv_u, degree, degree, degree, control_points, { elems_to_refine } );
            return drcd;
        }
        case TestCase::Case2p1KittyCornersIntersection_2d:
        {
            const size_t degree = 4;
            const size_t n_elems_1d = 2 * degree + 2;
            const basis::KnotVector kv_s = basis::unitIntervalKnotVectorWithNElems( n_elems_1d, degree );

            std::cout << n_elems_1d*grevillePoints( kv_s, degree ).transpose() << std::endl;
            const basis::KnotVector kv_t = kv_s;
            const Eigen::MatrixXd control_points = util::tensorProduct( { grevillePoints( kv_s, degree ), grevillePoints( kv_t, degree ) } ).transpose();

            std::vector<std::pair<size_t, size_t>> elems_to_refine;
            for( size_t i = 0; i <= degree; i++ )
            {
                for( size_t j = 0; j <= 2*degree; j++ )
                {
                    elems_to_refine.push_back( { j, i } );
                    if( j <= degree )
                        elems_to_refine.push_back( { n_elems_1d - 1 - j, n_elems_1d - 1 - i } );
                }
            }

            DeRhamComplexHierarchicalDiscretization drcd( kv_s, kv_t, degree, degree, control_points, { elems_to_refine } );
            return drcd;
        }
    }
}