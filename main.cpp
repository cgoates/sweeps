#include <iostream>
#include "SimplicialComplex.hpp"

int main()
{
    std::cout << "Hello World" << std::endl;

    std::vector<Simplex> simplices;
    simplices.push_back({1, 2, 3, 9});
    simplices.push_back({2, 4, 3, 9});
    simplices.push_back({1, 6, 2, 9});
    simplices.push_back({1, 5, 6, 9});
    simplices.push_back({3, 4, 8, 9});
    simplices.push_back({7, 3, 8, 9});
    simplices.push_back({2, 6, 4, 9});
    simplices.push_back({6, 8, 4, 9});
    simplices.push_back({3, 7, 5, 9});
    simplices.push_back({1, 3, 5, 9});
    simplices.push_back({6, 5, 8, 9});
    simplices.push_back({8, 5, 7, 9});

    const SimplicialComplex tets( simplices );

    std::cout << simplices.size() << std::endl;
}