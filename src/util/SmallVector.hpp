#pragma once
#include <ecpp/static_vector.hpp>

template<typename T, size_t MAX_SIZE>
using SmallVector = ecpp::static_vector<T, MAX_SIZE>;