/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2018-2023 OpenCFD Ltd.
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

#include "pointBoundaryMesh.H"
#include "polyBoundaryMesh.H"
#include "facePointPatch.H"
#include "pointMesh.H"
#include "PstreamBuffers.H"
#include "lduSchedule.H"
#include "globalMeshData.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::pointBoundaryMesh::addPatches(const polyBoundaryMesh& pbm)
{
    // Set boundary patches
    pointPatchList& patches = *this;

    patches.resize_null(pbm.size());

    forAll(patches, patchi)
    {
        // NB: needs ptr() to get *pointPatch instead of *facePointPatch
        patches.set(patchi, facePointPatch::New(pbm[patchi], *this).ptr());
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pointBoundaryMesh::pointBoundaryMesh
(
    const pointMesh& m,
    const polyBoundaryMesh& pbm
)
:
    pointPatchList(),
    mesh_(m)
{
    addPatches(pbm);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelList Foam::pointBoundaryMesh::indices
(
    const wordRe& matcher,
    const bool useGroups
) const
{
    return mesh().boundaryMesh().indices(matcher, useGroups);
}


Foam::labelList Foam::pointBoundaryMesh::indices
(
    const wordRes& matcher,
    const bool useGroups
) const
{
    return mesh().boundaryMesh().indices(matcher, useGroups);
}


Foam::labelList Foam::pointBoundaryMesh::indices
(
    const wordRes& select,
    const wordRes& ignore,
    const bool useGroups
) const
{
    return mesh().boundaryMesh().indices(select, ignore, useGroups);
}


Foam::label Foam::pointBoundaryMesh::findPatchID(const word& patchName) const
{
    return mesh().boundaryMesh().findPatchID(patchName);
}


void Foam::pointBoundaryMesh::calcGeometry()
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        pBufs.commsType() == Pstream::commsTypes::blocking
     || pBufs.commsType() == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initGeometry(pBufs);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).calcGeometry(pBufs);
        }
    }
    else if (pBufs.commsType() == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;

            if (schedEval.init)
            {
                operator[](patchi).initGeometry(pBufs);
            }
            else
            {
                operator[](patchi).calcGeometry(pBufs);
            }
        }
    }
}


void Foam::pointBoundaryMesh::movePoints(const pointField& p)
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        pBufs.commsType() == Pstream::commsTypes::blocking
     || pBufs.commsType() == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initMovePoints(pBufs, p);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).movePoints(pBufs, p);
        }
    }
    else if (pBufs.commsType() == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;

            if (schedEval.init)
            {
                operator[](patchi).initMovePoints(pBufs, p);
            }
            else
            {
                operator[](patchi).movePoints(pBufs, p);
            }
        }
    }
}


void Foam::pointBoundaryMesh::updateMesh()
{
    PstreamBuffers pBufs(Pstream::defaultCommsType);

    if
    (
        pBufs.commsType() == Pstream::commsTypes::blocking
     || pBufs.commsType() == Pstream::commsTypes::nonBlocking
    )
    {
        forAll(*this, patchi)
        {
            operator[](patchi).initUpdateMesh(pBufs);
        }

        pBufs.finishedSends();

        forAll(*this, patchi)
        {
            operator[](patchi).updateMesh(pBufs);
        }
    }
    else if (pBufs.commsType() == Pstream::commsTypes::scheduled)
    {
        const lduSchedule& patchSchedule = mesh().globalData().patchSchedule();

        // Dummy.
        pBufs.finishedSends();

        for (const auto& schedEval : patchSchedule)
        {
            const label patchi = schedEval.patch;

            if (schedEval.init)
            {
                operator[](patchi).initUpdateMesh(pBufs);
            }
            else
            {
                operator[](patchi).updateMesh(pBufs);
            }
        }
    }
}


// ************************************************************************* //
