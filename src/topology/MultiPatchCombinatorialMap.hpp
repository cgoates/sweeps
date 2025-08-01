#pragma once
#include <CombinatorialMap.hpp>
#include <TPCombinatorialMap.hpp>
#include <map>
#include <DartRange.hpp>

namespace topology
{
    class MultiPatchCombinatorialMap : public CombinatorialMap
    {
        public:
        struct ConstituentSide
        {
            size_t constituent_id;
            size_t side_id;

            bool operator<( const ConstituentSide& o ) const;
        };

        enum class TPPermutation : Dart::IndexType
        {
            ZeroToZero,
            ZeroToOne,
            ZeroToTwo,
            ZeroToThree,
            Flip1d
        };

        using InternalConnectionsMap =
            std::map<ConstituentSide, std::pair<TPPermutation, ConstituentSide>>;

        MultiPatchCombinatorialMap( const std::vector<std::shared_ptr<const TPCombinatorialMap>>& constituents,
                                    const std::map<std::pair<size_t, Dart>, std::pair<size_t, Dart>>& connections );
        MultiPatchCombinatorialMap( const std::vector<std::shared_ptr<const TPCombinatorialMap>>& constituents,
                                    const InternalConnectionsMap& connections );

        virtual ~MultiPatchCombinatorialMap() = default;

        virtual std::optional<Dart> phi( const int i, const Dart& d ) const override;

        virtual Dart::IndexType maxDartId() const override;

        virtual uint dim() const override;

        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override;

        virtual bool iterateCellsWhile( const uint cell_dim,
                                        const std::function<bool( const Cell& )>& callback ) const override;

        virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override;

        virtual std::optional<size_t> cellCount( const uint ) const override;

        std::pair<size_t, Dart> toLocalDart( const Dart& global_dart ) const;
        Dart toGlobalDart( const size_t patch_id, const Dart& local_dart ) const;

        const std::vector<std::shared_ptr<const TPCombinatorialMap>>& constituents() const { return mSubMaps; }
        const InternalConnectionsMap& connections() const { return mInterMapConnections; }

        private:
        std::vector<std::shared_ptr<const TPCombinatorialMap>> mSubMaps;
        std::map<ConstituentSide, std::pair<TPPermutation, ConstituentSide>> mInterMapConnections;
        const DartRanges mRanges;
    };

    std::ostream& operator<<( std::ostream& o, const MultiPatchCombinatorialMap::ConstituentSide& cs );
    std::ostream& operator<<( std::ostream& o, const MultiPatchCombinatorialMap::TPPermutation& cs );

    /// @brief Computes the connections of a 2d multipatch cmap swept along the third dimension.
    /// That is, if you tensor product each constituent of the 2d multipatch with the same 1d line
    /// and then want to connect the resulting 2d TP cmaps consistent with the 2d multipatch, this
    /// function will compute the necessary connections.
    /// @param connections_2d The connections from the 2d multipatch cmap.
    /// @return The connections for the 3d multipatch cmap.
    MultiPatchCombinatorialMap::InternalConnectionsMap
        connectionsOfSweptMultipatch( const MultiPatchCombinatorialMap::InternalConnectionsMap& connections_2d );

    /// @brief Computes the hex or quad layout of a multipatch combinatorial map by replacing each patch with a single
    /// element.
    MultiPatchCombinatorialMap blockLayout( const MultiPatchCombinatorialMap& cmap );
} // namespace topology