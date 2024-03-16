#include <array>
#include <span>

template<typename T, size_t MAX_SIZE>
class SmallVector
{
    public:
    SmallVector( const std::span<T>& s ) :
        mSize( s.size() )
    {
        for( size_t i = 0; i < s.size(); i++ ) push_back( s.at( i ) );
    }

    SmallVector() : mSize( 0 ) {}

    const T& at( const size_t n ) const
    {
        if( n >= mSize ) throw std::out_of_range( "Input index is greater than end of SmallVector" );
        return mData.at( n );
    }

    T& at( const size_t n )
    {
        if( n >= mSize ) throw std::out_of_range( "Input index is greater than end of SmallVector" );
        return mData.at( n );
    }

    void push_back( const T& t )
    {
        if( mSize == MAX_SIZE )
        {
            for( size_t i = 0; i < mSize; i++ ) std::cout << at( i ) << ", ";
            std::cout << std::endl;
            throw std::out_of_range( "push_back into full SmallVector" );
        }
        mData.at( mSize++ ) = t;
    }

    size_t size() const { return mSize; }

    std::array<T,MAX_SIZE>::iterator begin()
    {
        return mData.begin();
    }

    std::array<T,MAX_SIZE>::iterator end()
    {
        return std::next( mData.begin(), mSize );
    }

    std::array<T,MAX_SIZE>::const_iterator begin() const
    {
        return mData.begin();
    }

    std::array<T,MAX_SIZE>::const_iterator end() const
    {
        return std::next( mData.begin(), mSize );
    }

    bool operator==( const SmallVector<T, MAX_SIZE>& o ) const
    {
        if( size() != o.size() ) return false;
        for( size_t i = 0; i <= o.dim(); i++ )
        {
            if( at( i ) != o.at( i ) ) return false;
        }
        return true;
    }

    private:
    // TODO: Add optional std::vector for when it overflows
    size_t mSize;
    std::array<T, MAX_SIZE> mData;
};