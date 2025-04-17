#include <SmallVector.hpp>
#include <vector>

template <typename T, std::size_t N> class GrowableVector
{
    private:
    enum class StorageType
    {
        Static,
        Dynamic
    };

    StorageType mActiveStorage = StorageType::Static;
    SmallVector<T, N> mStaticStorage;
    std::vector<T> mDynamicStorage;

    // Convert from static to dynamic storage
    void switchToDynamic()
    {
        if( mActiveStorage == StorageType::Static )
        {
            mDynamicStorage.reserve( mStaticStorage.size() * 2 );
            for( size_t i = 0; i < mStaticStorage.size(); ++i )
            {
                mDynamicStorage.push_back( std::move( mStaticStorage.at(i) ) );
            }
            mStaticStorage.clear();
            mActiveStorage = StorageType::Dynamic;
        }
    }

    public:
    // Constructor
    GrowableVector() = default;

    // Constructor with initial capacity
    explicit GrowableVector( std::size_t capacity )
    {
        if( capacity > N )
        {
            mActiveStorage = StorageType::Dynamic;
            mDynamicStorage.reserve( capacity );
        }
    }

    // Copy constructor
    GrowableVector( const GrowableVector& other ) : mActiveStorage( other.mActiveStorage )
    {
        if( other.mActiveStorage == StorageType::Static )
        {
            mStaticStorage = other.mStaticStorage;
        }
        else
        {
            mDynamicStorage = other.mDynamicStorage;
        }
    }

    // Copy assignment
    GrowableVector& operator=( const GrowableVector& other )
    {
        if( this != &other )
        {
            if( other.mActiveStorage == StorageType::Static )
            {
                mStaticStorage = other.mStaticStorage;
                mDynamicStorage.clear();
                mActiveStorage = StorageType::Static;
            }
            else
            {
                mDynamicStorage = other.mDynamicStorage;
                mStaticStorage.clear();
                mActiveStorage = StorageType::Dynamic;
            }
        }
        return *this;
    }

    T& at( std::size_t index )
    {
        if( mActiveStorage == StorageType::Static )
        {
            return mStaticStorage.at( index );
        }
        else
        {
            return mDynamicStorage.at( index );
        }
    }

    const T& at( std::size_t index ) const
    {
        if( mActiveStorage == StorageType::Static )
        {
            return mStaticStorage.at( index );
        }
        else
        {
            return mDynamicStorage.at( index );
        }
    }

    T& front()
    {
        if( mActiveStorage == StorageType::Static )
        {
            return mStaticStorage.front();
        }
        else
        {
            return mDynamicStorage.front();
        }
    }

    const T& front() const
    {
        if( mActiveStorage == StorageType::Static )
        {
            return mStaticStorage.front();
        }
        else
        {
            return mDynamicStorage.front();
        }
    }

    T& back()
    {
        if( mActiveStorage == StorageType::Static )
        {
            return mStaticStorage.back();
        }
        else
        {
            return mDynamicStorage.back();
        }
    }

    const T& back() const
    {
        if( mActiveStorage == StorageType::Static )
        {
            return mStaticStorage.back();
        }
        else
        {
            return mDynamicStorage.back();
        }
    }

    // Capacity
    std::size_t size() const
    {
        if( mActiveStorage == StorageType::Static )
        {
            return mStaticStorage.size();
        }
        else
        {
            return mDynamicStorage.size();
        }
    }

    std::size_t capacity() const
    {
        if( mActiveStorage == StorageType::Static )
        {
            return N;
        }
        else
        {
            return mDynamicStorage.capacity();
        }
    }

    bool empty() const
    {
        if( mActiveStorage == StorageType::Static )
        {
            return mStaticStorage.empty();
        }
        else
        {
            return mDynamicStorage.empty();
        }
    }

    void reserve( std::size_t new_capacity )
    {
        if( new_capacity > N && mActiveStorage == StorageType::Static )
        {
            switchToDynamic();
            mDynamicStorage.reserve( new_capacity );
        }
        else if( mActiveStorage == StorageType::Dynamic )
        {
            mDynamicStorage.reserve( new_capacity );
        }
    }

    // Modifiers
    void push_back( const T& value )
    {
        if( mActiveStorage == StorageType::Static )
        {
            if( mStaticStorage.size() < N )
            {
                mStaticStorage.push_back( value );
            }
            else
            {
                switchToDynamic();
                mDynamicStorage.push_back( value );
            }
        }
        else
        {
            mDynamicStorage.push_back( value );
        }
    }

    void push_back( T&& value )
    {
        if( mActiveStorage == StorageType::Static )
        {
            if( mStaticStorage.size() < N )
            {
                mStaticStorage.push_back( std::move( value ) );
            }
            else
            {
                switchToDynamic();
                mDynamicStorage.push_back( std::move( value ) );
            }
        }
        else
        {
            mDynamicStorage.push_back( std::move( value ) );
        }
    }

    template <typename... Args> T& emplace_back( Args&&... args )
    {
        if( mActiveStorage == StorageType::Static )
        {
            if( mStaticStorage.size() < N )
            {
                return mStaticStorage.emplace_back( std::forward<Args>( args )... );
            }
            else
            {
                switchToDynamic();
                return mDynamicStorage.emplace_back( std::forward<Args>( args )... );
            }
        }
        else
        {
            return mDynamicStorage.emplace_back( std::forward<Args>( args )... );
        }
    }

    void pop_back()
    {
        if( mActiveStorage == StorageType::Static )
        {
            mStaticStorage.pop_back();
        }
        else
        {
            mDynamicStorage.pop_back();
        }
    }

    void clear()
    {
        if( mActiveStorage == StorageType::Static )
        {
            mStaticStorage.clear();
        }
        else
        {
            mDynamicStorage.clear();
        }
    }

    // Iterators
    T* begin()
    {
        if( mActiveStorage == StorageType::Static )
        {
            return mStaticStorage.begin();
        }
        else
        {
            return mDynamicStorage.data();
        }
    }

    const T* begin() const
    {
        if( mActiveStorage == StorageType::Static )
        {
            return mStaticStorage.begin();
        }
        else
        {
            return mDynamicStorage.data();
        }
    }

    T* end()
    {
        if( mActiveStorage == StorageType::Static )
        {
            return mStaticStorage.end();
        }
        else
        {
            return begin() + mDynamicStorage.size();
        }
    }

    const T* end() const
    {
        if( mActiveStorage == StorageType::Static )
        {
            return mStaticStorage.end();
        }
        else
        {
            return begin() + mDynamicStorage.size();
        }
    }

    // Shrink to fit - move back to static if possible
    void shrink_to_fit()
    {
        if( mActiveStorage == StorageType::Dynamic && mDynamicStorage.size() <= N )
        {
            SmallVector<T, N> temp_static;
            for( const auto& item : mDynamicStorage )
            {
                temp_static.push_back( item );
            }
            mDynamicStorage.clear();
            mStaticStorage = std::move( temp_static );
            mActiveStorage = StorageType::Static;
        }
        else if( mActiveStorage == StorageType::Dynamic )
        {
            mDynamicStorage.shrink_to_fit();
        }
    }

    // Check if we're using dynamic storage
    bool is_dynamic() const { return mActiveStorage == StorageType::Dynamic; }
};