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
#include "momentumWeightedReconstruct.H"
#include "fvcGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::momentumWeightedReconstruct::correction
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
                "momentumWeightedReconstruct::correction(" + vf.name() + ')',
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

    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    const volVectorField& C = mesh_.C();
    const surfaceVectorField& Cf = mesh_.Cf();

    IStringStream scheme("VKLimited leastSquares 1.0");

    tmp<volVectorField> tgradVf
    (
        Foam::fv::gradScheme<scalar>::New(mesh_, scheme)().grad(vf)
    );

    volVectorField& gradVf = tgradVf.ref();

    // Adjust pressure gradient at porous cells
    /*if (hasActivePorousZone_)
    {
        const volVectorField& porousGradp
        (
            mesh_.lookupObject<volVectorField>("porousGradp")
        );

        const labelList& porousCells = fvOptions_.activePorousCells();

        forAll(porousCells, i)
        {
            const label& celli(porousCells[i]);

            gradVf[celli] = porousGradp[celli];
        }
    }

    gradVf.correctBoundaryConditions();*/

    const scalarField& w = weights_.internalField();

    if (hasActivePorousZone_)
    {
        const surfaceScalarField& pFs(fvOptions_.porousFaces());

        forAll(w, facei)
        {
            if (!pFs[facei])
            {
                const label& own = owner[facei];
                const label& nei = neighbour[facei];

                sfCorr[facei] =
                    w[facei]*((Cf[facei] - C[own]) & gradVf[own])
                  + (1 - w[facei])*((Cf[facei] - C[nei]) & gradVf[nei]);
            } 
        }
    }
    else
    {
        forAll(w, facei)
        {
            const label& own = owner[facei];
            const label& nei = neighbour[facei];

            sfCorr[facei] =
                w[facei]*((Cf[facei] - C[own]) & gradVf[own])
              + (1 - w[facei])*((Cf[facei] - C[nei]) & gradVf[nei]);
        }
    }

    typename surfaceScalarField::Boundary& bSfCorr = sfCorr.boundaryFieldRef();

    forAll(bSfCorr, patchi)
    {
        fvsPatchField<scalar>& pSfCorr = bSfCorr[patchi];

        const polyPatch& p = mesh_.boundary()[patchi].patch();

        const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();

        const vectorField& pCf = Cf.boundaryField()[patchi];

        const scalarField& pWeights = weights_.boundaryField()[patchi];

        if (p.coupled())
        {
            const vectorField pGradVfNei
            (
                gradVf.boundaryField()[patchi].patchNeighbourField()
            );

            // Build the d-vectors
            vectorField pd(Cf.boundaryField()[patchi].patch().delta());

            if (hasActivePorousZone_)
            {            
                auto& pPorousFaces = 
                    fvOptions_.porousFaces().boundaryField()[patchi];

                forAll(pOwner, facei)
                {
                    if (pPorousFaces[facei])
                    {
                        label own = pOwner[facei];

                        pSfCorr[facei] =
                            pWeights[facei]
                            *((pCf[facei] - C[own]) & gradVf[own])
                          + (1 - pWeights[facei])
                            *(
                                 (pCf[facei] - pd[facei] - C[own]) 
                               & pGradVfNei[facei]
                             );
                    }
                }
            }
            else
            {
                forAll(pOwner, facei)
                {
                    label own = pOwner[facei];

                    pSfCorr[facei] =
                        pWeights[facei]*((pCf[facei] - C[own]) & gradVf[own])
                      + (1 - pWeights[facei])
                        *(
                             (pCf[facei] - pd[facei] - C[own]) 
                           & pGradVfNei[facei]
                         );
                }
            }
        }
        /*else if(!isA<emptyPolyPatch>(p) && hasActivePorousZone_)
        {
            const scalarField& pPorousFaces = 
                fvOptions_.porousFaces().boundaryField()[patchi];

            forAll(pOwner, facei)
            {
                if (pPorousFaces[facei])
                {
                    label own = pOwner[facei];

                    pSfCorr[facei] = ((pCf[facei] - C[own]) & gradVf[own]);
                }
            }
        }*/
        //else if (!vf.boundaryField()[patchi].fixesValue())
        //{
        //    // For patches that do not fix the value, calculate
        //    // extrapolated field
        //    forAll(pOwner, facei)
        //    {
        //        label own = pOwner[facei];
        //
        //        vector df = pCf[facei] - C[own];
        //
        //        pSfCorr[facei] = vf[own] - vf.boundaryField()[patchi][facei]
        //            + pWeights[facei]*(df & gradVf[own]);
        //    }
        //}
    }

    return tsfCorr;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePressureInterpolationScheme(momentumWeightedReconstruct)
}

// ************************************************************************* //
