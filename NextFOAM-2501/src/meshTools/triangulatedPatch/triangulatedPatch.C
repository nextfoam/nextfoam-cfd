/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "triangulatedPatch.H"
#include "triPointRef.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::triangulatedPatch::randomPoint
(
    Random& rnd,
    const scalar c,
    point& result,
    label& facei,
    label& celli
) const
{
    result = point::min;
    facei = -1;
    celli = -1;

    if (triWght_.empty() || c < triWght_.front() || c > triWght_.back())
    {
        return false;
    }

    // Find corresponding decomposed face triangle
    // Note: triWght_ is sized nTri+1 (zero added at start)
    label trii = 0;
    for (label i = 0; i < triWght_.size() - 1; ++i)
    {
        if (c > triWght_[i] && c <= triWght_[i+1])
        {
            trii = i;
            break;
        }
    }

    // Find random point in triangle
    const pointField& points = patch_.points();
    const face& tf = triFace_[trii];
    const triPointRef tri(points[tf[0]], points[tf[1]], points[tf[2]]);

    result = tri.randomPoint(rnd);
    facei = triToFace_[trii];
    celli = patch_.faceCells()[facei];

    if (perturbTol_ > 0)
    {
        const polyMesh& mesh = patch_.boundaryMesh().mesh();
        const point& cc = mesh.cellCentres()[celli];

        // Normal points out of domain => subtract correction
        const vector& n = patch_.faceNormals()[facei];
        result -= perturbTol_*n*mag(n & (cc - result));

        // Reset facei - point no longer resides on the face
        facei = -1;
    }

    return true;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triangulatedPatch::triangulatedPatch
(
    const polyPatch& patch,
    const scalar perturbTol
)
:
    patch_(patch),
    perturbTol_(perturbTol),
    triFace_(),
    triToFace_(),
    triWght_()
{
    update();
}


Foam::triangulatedPatch::triangulatedPatch
(
    const polyMesh& mesh,
    const word& patchName,
    const scalar perturbTol
)
:
    triangulatedPatch(mesh.boundaryMesh()[patchName], perturbTol)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::triangulatedPatch::update()
{
    const pointField& points = patch_.points();

    // Triangulate the patch faces and create addressing
    DynamicList<label> triToFace(2*patch_.size());
    DynamicList<face> triFace(2*patch_.size());
    DynamicList<scalar> triWght(2*patch_.size());
    DynamicList<face> tris(8);

    // Set zero value at the start of the tri area/weight list
    triWght.push_back(0);

    forAll(patch_, facei)
    {
        const face& f = patch_[facei];

        tris.clear();
        f.triangles(points, tris);

        for (const auto& t : tris)
        {
            triToFace.push_back(facei);
            triFace.push_back(t);
            triWght.push_back(t.mag(points));
        }
    }

    scalarList procSumWght(Pstream::nProcs()+1, Zero);
    procSumWght[Pstream::myProcNo()+1] = sum(triWght);
    Pstream::listCombineReduce(procSumWght, maxEqOp<scalar>());

    for (label i = 1; i < procSumWght.size(); ++i)
    {
        // Convert to cumulative
        procSumWght[i] += procSumWght[i-1];
    }

    const scalar offset = procSumWght[Pstream::myProcNo()];
    forAll(triWght, i)
    {
        if (i)
        {
            // Convert to cumulative
            triWght[i] += triWght[i-1];
        }

        // Apply processor offset
        triWght[i] += offset;
    }

    // Normalise
    const scalar sumWght = procSumWght.back();
    for (scalar& w : triWght)
    {
        w /= sumWght;
    }

    // Transfer to persistent storage
    triFace_.transfer(triFace);
    triToFace_.transfer(triToFace);
    triWght_.transfer(triWght);
}


bool Foam::triangulatedPatch::randomLocalPoint
(
    Random& rnd,
    point& result,
    label& facei,
    label& celli
) const
{
    const scalar c = rnd.position<scalar>(triWght_.front(), triWght_.back());

    return randomPoint(rnd, c, result, facei, celli);
}


bool Foam::triangulatedPatch::randomGlobalPoint
(
    Random& rnd,
    point& result,
    label& facei,
    label& celli
) const
{
    boolList valid(UPstream::nProcs(), false);
    valid[UPstream::myProcNo()] = randomLocalPoint(rnd, result, facei, celli);
    UPstream::listGatherValues(valid);

    forAll(valid, proci)
    {
        // Choose first valid processor
        if (valid[proci])
        {
            return (proci == UPstream::myProcNo());
        }
    }

    return false;
}


// ************************************************************************* //
