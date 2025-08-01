#pragma once
#include <CombinatorialMap.hpp>
#include <CombinatorialMap1d.hpp>
#include <SmallVector.hpp>

namespace topology
{
    /// @brief A combinatorial map formed from the tensor product of another topology (the 'source') and a
    /// one-dimensional topology (the 'line') to create a ( source.dim() + 1 )-dimensional topology.
    ///
    /// If source.dim() is 1, each line dart/source dart combo becomes four darts.
    /// DP# means the dart has TPDartPos::DartPos#.  This is how we differentiate
    /// between the four darts coming from the same source and line dart.
    /// *  ------------*
    /// |      DP2
    /// |              |
    /// | DP3      DP1 |
    /// |              |
    ///        DP0     |
    /// *------------  *
    ///
    /// If source.dim() is 2, each line dart/source dart combo becomes 6 darts.
    /// Two of these are on faces perpendicular to the sweep, and four form a
    /// face that is parallel to the sweep.
    ///
    ///        DP5
    ///   ------------*
    ///  *------------
    ///        DP3     *
    /// |              |
    /// | DP2      DP4 |
    /// |              |
    /// *      DP1     |
    ///   ------------*
    ///  *------------
    ///        DP0
    ///
    class TPCombinatorialMap : public CombinatorialMap
    {
        public:
        TPCombinatorialMap( const std::shared_ptr<const CombinatorialMap>& source,
                            const std::shared_ptr<const CombinatorialMap1d>& line );
        virtual ~TPCombinatorialMap() = default;

        virtual std::optional<Dart> phi( const int i, const Dart& d ) const override;

        virtual Dart::IndexType maxDartId() const override;

        virtual uint dim() const override;

        virtual bool iterateDartsWhile( const std::function<bool( const Dart& )>& callback ) const override;

        virtual bool iterateCellsWhile( const uint cell_dim,
                                        const std::function<bool( const Cell& )>& callback ) const override;

        virtual std::optional<IndexingFunc> indexing( const uint cell_dim ) const override;

        virtual std::optional<size_t> cellCount( const uint ) const override;

        const CombinatorialMap& sourceCMap() const { return *mSource; }
        const CombinatorialMap1d& lineCMap() const { return *mLine; }
        const std::shared_ptr<const CombinatorialMap>& sourceCMapPtr() const { return mSource; }
        const std::shared_ptr<const CombinatorialMap1d>& lineCMapPtr() const { return mLine; }

        enum class TPDartPos : Dart::IndexType
        {
            DartPos0,
            DartPos1,
            DartPos2,
            DartPos3,
            DartPos4,
            DartPos5
        };

        Dart flatten( const Dart& source_dart, const Dart& line_dart, const TPDartPos& pos ) const;
        std::tuple<Dart, Dart, TPDartPos> unflatten( const Dart& d ) const;

        private:

        Dart::IndexType dartsPerSourceDart() const { return 2 * dim(); }

        const std::shared_ptr<const CombinatorialMap> mSource;
        const std::shared_ptr<const CombinatorialMap1d> mLine;
    };

    SmallVector<std::shared_ptr<const CombinatorialMap1d>, 3> tensorProductComponentCMaps( const TPCombinatorialMap& tp_map );
    TPCombinatorialMap tensorProductCMapFromComponents( const SmallVector<std::shared_ptr<const CombinatorialMap1d>, 3>& components );

    struct FullyUnflattenedDart
    {
        FullyUnflattenedDart( const SmallVector<Dart, 3>& darts,
                              const SmallVector<TPCombinatorialMap::TPDartPos, 2>& pos )
            : unflat_darts( darts ), dart_pos( pos )
        {}
        FullyUnflattenedDart( const SmallVector<Dart, 3>& darts )
            : FullyUnflattenedDart( darts,
                                    SmallVector<TPCombinatorialMap::TPDartPos, 2>(
                                        darts.size() - 1, TPCombinatorialMap::TPDartPos::DartPos0 ) )
        {}
        FullyUnflattenedDart() = default;
        SmallVector<Dart, 3> unflat_darts;
        SmallVector<TPCombinatorialMap::TPDartPos, 2> dart_pos;
    };

    FullyUnflattenedDart unflattenFull( const TPCombinatorialMap& cmap, const Dart& d );
    Dart flattenFull( const TPCombinatorialMap& cmap, const FullyUnflattenedDart& unflat_d );

    using CellOrEndVertex = std::optional<Cell>;
    /// @brief Returns the two cells that the given TP cell is a tensor product of.
    /// Careful in using this, doesn't care about dart orientation or which element the darts end up in.
    /// @param cmap The TP CMap.
    /// @param c    The TP Cell to unflatten.
    /// @return     A pair { source_cell, line_cell }. source_cell.dim() + line_cell.dim() = c.dim()
    std::pair<CellOrEndVertex, CellOrEndVertex> unflattenCell( const TPCombinatorialMap& cmap, const Cell& c );
}; // namespace topology