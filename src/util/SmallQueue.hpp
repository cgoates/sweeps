#include <SmallVector.hpp>
#include <queue>
#include <deque>

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

template <typename T, size_t N> class GrowableQueue
{
    public:
    GrowableQueue() : mUsingSmallVector( true ) {}

    void push( const T& t )
    {
        if( mUsingSmallVector )
        {
            if( mSmallData.size() < N )
            {
                mSmallData.push_back( t );
            }
            else
            {
                switchToQueue();
                mQueue.push( t );
            }
        }
        else
            mQueue.push( t );
    }

    T pop()
    {
        if( mUsingSmallVector )
        {
            return mSmallData.at( mFrontPos++ );
        }
        else
        {
            T value = mQueue.front();
            mQueue.pop();
            return value;
        }
    }

    bool empty() const
    {
        if( mUsingSmallVector )
        {
            return mFrontPos == mSmallData.size();
        }
        else
        {
            return mQueue.empty();
        }
    }

    private:
    void switchToQueue()
    {
        for( size_t i = mFrontPos; i < mSmallData.size(); ++i )
        {
            mQueue.push( mSmallData.at( i ) );
        }

        mSmallData.clear();
        mFrontPos = 0;

        mUsingSmallVector = false;
    }

    bool mUsingSmallVector;
    size_t mFrontPos = 0;
    SmallVector<T, N> mSmallData;
    std::queue<T, std::deque<T>> mQueue;
};