/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "fvMesh.H"
#include "momentumWeighted.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::momentumWeighted::correction
(
     const volScalarField& vf
) const
{
    tmp<surfaceScalarField> tsfCorr
    (
        new surfaceScalarField
        (
            IOobject
            (
                "momentumWeighted::correction(" + vf.name() + ')',
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            mesh_,
            dimensionedScalar(vf.name(), vf.dimensions(), Zero)
        )
    );

    surfaceScalarField& sfCorr = tsfCorr.ref();

    const surfaceScalarField& pFs = fvOptions_.porousFaces();

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const volVectorField& C = mesh_.C();
    const surfaceVectorField& Cf = mesh_.Cf();

    const volVectorField& gradVf
    (
        mesh_.lookupObject<volVectorField>("porousGradp")
    );

    const scalarField& w = weights_.internalField();
 
    forAll(pFs, facei)
    {
        if (pFs[facei])
        {
            const label& own = owner[facei];
            const label& nei = neighbour[facei];

            sfCorr[facei] =
            w[facei]
            *(
                ((Cf[facei] - C[own]) & gradVf[own])
              + ((Cf[facei] - C[nei]) & gradVf[nei])
             );
        }
    }

    typename surfaceScalarField::Boundary& bSfCorr = sfCorr.boundaryFieldRef();

    forAll(bSfCorr, patchi)
    {
        fvsPatchField<scalar>& pSfCorr = bSfCorr[patchi];

        const polyPatch& p = mesh_.boundary()[patchi].patch();

        const labelUList& pOwner = p.faceCells();

        const vectorField& pCf = Cf.boundaryField()[patchi];

        const scalarField& pPorousFaces = pFs.boundaryField()[patchi];

        const scalarField& pWeights = weights_.boundaryField()[patchi];

        if (p.coupled())
        {
            const vectorField pGradVfNei
            (
                gradVf.boundaryField()[patchi].patchNeighbourField()
            );

            // Build the d-vectors
            vectorField pd(Cf.boundaryField()[patchi].patch().delta());

            forAll(pOwner, facei)
            {
                if (pPorousFaces[facei])
                {
                    label own = pOwner[facei];
                    vector Cnei = C[own] + pd[facei];

                    pSfCorr[facei] =
                        pWeights[facei]
                        *(
                            ((pCf[facei] - C[own]) & gradVf[own])
                          + ((pCf[facei] - Cnei) & pGradVfNei[facei])
                         );
                }
            }
        }
        else if(!isA<emptyPolyPatch>(p))
        {
            forAll(pOwner, facei)
            {
                if (pPorousFaces[facei])
                {
                    label own = pOwner[facei];

                    pSfCorr[facei] = ((pCf[facei] - C[own]) & gradVf[own]);
                }
            }
        }
    }

    return tsfCorr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePressureInterpolationScheme(momentumWeighted)
}

// ************************************************************************* //
