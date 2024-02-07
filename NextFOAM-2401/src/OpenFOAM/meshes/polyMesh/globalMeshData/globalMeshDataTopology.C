/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
    Copyright (C) 2015-2023 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "globalMeshData.H"
#include "globalIndex.H"
#include "CompactListList.H"
#include "processorPolyPatch.H"
#include "syncTools.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

void Foam::globalMeshData::calcCellCells
(
    const polyMesh& mesh,
    const labelList& agglom,
    const label nLocalCoarse,
    const bool parallel,
    CompactListList<label>& cellCells
)
{
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeigh = mesh.faceNeighbour();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // FUTURE? treat empty agglomeration like an identity map

    // Global cell numbers (agglomerated numbering)
    const label myProci = UPstream::myProcNo(UPstream::worldComm);
    const globalIndex globalAgglom(nLocalCoarse, UPstream::worldComm, parallel);


    // The agglomerated owner per boundary faces (global numbering)
    // from the other side (of coupled patches)

    labelList globalNeighbour(agglom, pbm.faceOwner());
    globalAgglom.inplaceToGlobal(myProci, globalNeighbour);
    syncTools::swapBoundaryFaceList(mesh, globalNeighbour);


    // Count number of faces (internal + coupled)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Number of faces per coarse cell
    labelList nFacesPerCell(nLocalCoarse, Zero);

    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const label own = agglom[faceOwner[facei]];
        const label nei = agglom[faceNeigh[facei]];

        ++nFacesPerCell[own];
        ++nFacesPerCell[nei];
    }

    for (const polyPatch& pp : pbm)
    {
        if (pp.coupled() && (parallel || !isA<processorPolyPatch>(pp)))
        {
            label bFacei = pp.start()-mesh.nInternalFaces();

            for (const label facei : pp.range())
            {
                const label own = agglom[faceOwner[facei]];
                const label globalNei = globalNeighbour[bFacei];

                if
                (
                   !globalAgglom.isLocal(myProci, globalNei)
                 || globalAgglom.toLocal(myProci, globalNei) != own
                )
                {
                    ++nFacesPerCell[own];
                }

                ++bFacei;
            }
        }
    }


    // Fill in offset and data
    // ~~~~~~~~~~~~~~~~~~~~~~~

    cellCells.resize_nocopy(nFacesPerCell);

    nFacesPerCell = 0;  // Restart the count

    labelList& m = cellCells.values();
    const labelList& offsets = cellCells.offsets();

    // For internal faces is just offsetted owner and neighbour
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const label own = agglom[faceOwner[facei]];
        const label nei = agglom[faceNeigh[facei]];

        const label ownIndex = offsets[own] + nFacesPerCell[own]++;
        const label neiIndex = offsets[nei] + nFacesPerCell[nei]++;

        m[ownIndex] = globalAgglom.toGlobal(myProci, nei);
        m[neiIndex] = globalAgglom.toGlobal(myProci, own);
    }

    // For boundary faces is offsetted coupled neighbour
    for (const polyPatch& pp : pbm)
    {
        if (pp.coupled() && (parallel || !isA<processorPolyPatch>(pp)))
        {
            label bFacei = pp.start()-mesh.nInternalFaces();

            for (const label facei : pp.range())
            {
                const label own = agglom[faceOwner[facei]];
                const label globalNei = globalNeighbour[bFacei];

                if
                (
                   !globalAgglom.isLocal(myProci, globalNei)
                 || globalAgglom.toLocal(myProci, globalNei) != own
                )
                {
                    const label ownIndex = offsets[own] + nFacesPerCell[own]++;

                    m[ownIndex] = globalNei;
                }

                ++bFacei;
            }
        }
    }


    // Check for duplicates connections between cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Done as postprocessing step since we now have cellCells.

    // NB: Because of agglomeration, self-connections will occur
    //     and must be filtered out.

    if (!cellCells.empty())
    {
        label newIndex = 0;
        labelHashSet nbrCells;

        labelList& m = cellCells.values();
        labelList& offsets = cellCells.offsets();

        label startIndex = offsets[0];

        const label nCellCells = cellCells.size();

        for (label celli = 0; celli < nCellCells; ++celli)
        {
            // Could be done as combination of std::sort, std::unique_copy
            // except need to block out 'self' and std::copy_if/remove_if
            // are undefined for overlapping regions

            const label self = globalAgglom.toGlobal(myProci, celli);

            nbrCells.clear();
            nbrCells.insert(self);

            const label endIndex = offsets[celli+1];

            for (label i = startIndex; i < endIndex; ++i)
            {
                if (nbrCells.insert(m[i]))
                {
                    m[newIndex] = m[i];
                    ++newIndex;
                }
            }
            startIndex = endIndex;
            offsets[celli+1] = newIndex;
        }

        m.resize(newIndex);
    }

    //forAll(cellCells, celli)
    //{
    //    Pout<< "Original: Coarse cell " << celli << endl;
    //    forAll(mesh.cellCells()[celli], i)
    //    {
    //        Pout<< "    nbr:" << mesh.cellCells()[celli][i] << endl;
    //    }
    //    Pout<< "Compacted: Coarse cell " << celli << endl;
    //    const labelUList cCells = cellCells[celli];
    //    forAll(cCells, i)
    //    {
    //        Pout<< "    nbr:" << cCells[i] << endl;
    //    }
    //}
}


void Foam::globalMeshData::calcCellCells
(
    const polyMesh& mesh,
    const labelList& agglom,
    const label nLocalCoarse,
    const bool parallel,
    CompactListList<label>& cellCells,
    CompactListList<scalar>& cellCellWeights
)
{
    const labelList& faceOwner = mesh.faceOwner();
    const labelList& faceNeigh = mesh.faceNeighbour();
    const polyBoundaryMesh& pbm = mesh.boundaryMesh();

    // FUTURE? treat empty agglomeration like an identity map

    // Global cell numbers (agglomerated numbering)
    const label myProci = UPstream::myProcNo(UPstream::worldComm);
    const globalIndex globalAgglom(nLocalCoarse, UPstream::worldComm, parallel);


    // The agglomerated owner per boundary faces (global numbering)
    // from the other side (of coupled patches)

    labelList globalNeighbour(agglom, pbm.faceOwner());
    globalAgglom.inplaceToGlobal(myProci, globalNeighbour);
    syncTools::swapBoundaryFaceList(mesh, globalNeighbour);


    // Count number of faces (internal + coupled)
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Number of faces per coarse cell
    labelList nFacesPerCell(nLocalCoarse, Zero);

    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const label own = agglom[faceOwner[facei]];
        const label nei = agglom[faceNeigh[facei]];

        ++nFacesPerCell[own];
        ++nFacesPerCell[nei];
    }

    for (const polyPatch& pp : pbm)
    {
        if (pp.coupled() && (parallel || !isA<processorPolyPatch>(pp)))
        {
            label bFacei = pp.start()-mesh.nInternalFaces();

            for (const label facei : pp.range())
            {
                const label own = agglom[faceOwner[facei]];
                const label globalNei = globalNeighbour[bFacei];

                if
                (
                   !globalAgglom.isLocal(myProci, globalNei)
                 || globalAgglom.toLocal(myProci, globalNei) != own
                )
                {
                    ++nFacesPerCell[own];
                }

                ++bFacei;
            }
        }
    }


    // Fill in offset and data
    // ~~~~~~~~~~~~~~~~~~~~~~~

    cellCells.resize_nocopy(nFacesPerCell);
    cellCellWeights.resize_nocopy(nFacesPerCell);

    nFacesPerCell = 0;  // Restart the count

    labelList& m = cellCells.values();
    scalarList& w = cellCellWeights.values();
    const labelList& offsets = cellCells.offsets();

    // For internal faces is just offsetted owner and neighbour
    for (label facei = 0; facei < mesh.nInternalFaces(); ++facei)
    {
        const label own = agglom[faceOwner[facei]];
        const label nei = agglom[faceNeigh[facei]];

        const label ownIndex = offsets[own] + nFacesPerCell[own]++;
        const label neiIndex = offsets[nei] + nFacesPerCell[nei]++;

        m[ownIndex] = globalAgglom.toGlobal(myProci, nei);
        m[neiIndex] = globalAgglom.toGlobal(myProci, own);

        w[ownIndex] = mag(mesh.faceAreas()[facei]);
        w[neiIndex] = w[ownIndex];
    }

    // For boundary faces is offsetted coupled neighbour
    for (const polyPatch& pp : pbm)
    {
        if (pp.coupled() && (parallel || !isA<processorPolyPatch>(pp)))
        {
            label bFacei = pp.start()-mesh.nInternalFaces();

            for (const label facei : pp.range())
            {
                const label own = agglom[faceOwner[facei]];
                const label globalNei = globalNeighbour[bFacei];

                if
                (
                   !globalAgglom.isLocal(myProci, globalNei)
                 || globalAgglom.toLocal(myProci, globalNei) != own
                )
                {
                    const label ownIndex = offsets[own] + nFacesPerCell[own]++;

                    m[ownIndex] = globalNei;
                    w[ownIndex] = mag(mesh.faceAreas()[facei]);
                }

                ++bFacei;
            }
        }
    }


    // Check for duplicates connections between cells
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Done as postprocessing step since we now have cellCells.

    // NB: Because of agglomeration, self-connections will occur
    //     and must be filtered out.

    if (!cellCells.empty())
    {
        label newIndex = 0;
        labelHashSet nbrCells;

        labelList& m = cellCells.values();
        scalarList& w = cellCellWeights.values();
        labelList& offsets = cellCells.offsets();

        label startIndex = offsets[0];

        const label nCellCells = cellCells.size();

        for (label celli = 0; celli < nCellCells; ++celli)
        {
            const label self = globalAgglom.toGlobal(myProci, celli);

            nbrCells.clear();
            nbrCells.insert(self);

            const label endIndex = offsets[celli+1];

            for (label i = startIndex; i < endIndex; ++i)
            {
                if (nbrCells.insert(m[i]))
                {
                    m[newIndex] = m[i];
                    w[newIndex] = w[i];
                    ++newIndex;
                }
            }
            startIndex = endIndex;
            offsets[celli+1] = newIndex;
        }

        m.resize(newIndex);
        w.resize(newIndex);

        // Weights has identical offsets as cellCells
        cellCellWeights.offsets() = cellCells.offsets();
    }
}


// ************************************************************************* //
