#include <SmallVector.hpp>

template<typename T, size_t MAX_SIZE>
class SmallQueue
{
    public:
    SmallQueue() {}

    void push( const T& t ) { mData.push_back( t ); }

    T pop() { return mData.at( front_pos++ ); }

    bool empty() const { return front_pos == mData.size(); }

    private:

    size_t front_pos = 0;
    SmallVector<T, MAX_SIZE> mData;
};