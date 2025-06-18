/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "spectralRadius.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
spectralRadius::spectralRadius
(
    const volScalarField& rho,
    const volVectorField& U,
    const fluidThermo& thermoPhysicalModel,
    const turbulenceModel& turbulenceModel
)
:
    runTime_(U.time()),
    mesh_(U.mesh()),

    rho_(rho),
    U_(U),

    p_(thermoPhysicalModel.p()),

    thermoPhysicalModel_(thermoPhysicalModel),
    turbulenceModel_(turbulenceModel),

    w_(1.1),
    spectralRadiiInv_
    (
        IOobject
        (
            "spectralRadiiInv",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(mag(U_))
    ),
    spectralRadiiVis_
    (
        IOobject
        (
            "spectralRadiiVis",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(mag(U_))
    ),
    spectralRadiiown_
    (
        IOobject
        (
            "spectralRadiiown",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(mag(U_))
    ),
    spectralRadiinei_
    (
        IOobject
        (
            "spectralRadiinei",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        linearInterpolate(mag(U_))
    )
{} 

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void spectralRadius::evaluateSpectralRadii
(
    scalar& spectralRadiiInv,
    scalar& spectralRadiiVis,
    scalar& spectralRadiiown,
    scalar& spectralRadiinei,

    vector& dLoc,
    const scalar alphaEffLeft,
    const scalar alphaEffRight,

    const scalar pLeft,
    const scalar pRight,
    const vector ULeft,
    const vector URight,
    const scalar rhoLeft,
    const scalar rhoRight,

    const scalar gammaLeft,
    const scalar gammaRight,

    const vector Sf,
    const scalar magSf
) const
{
    // normal vector
    vector normalVector = Sf/magSf;

    const scalar contrVLeft  = (ULeft & normalVector);
    const scalar contrVRight = (URight & normalVector);

    const scalar speedOfSoundLeft = sqrt(gammaLeft*pLeft/rhoLeft);
    const scalar speedOfSoundRight = sqrt(gammaRight*pRight/rhoRight);

    const scalar pAvg = 0.5*(pLeft + pRight);
    const scalar contrVavg = (0.5*(ULeft + URight)) & normalVector;
    const scalar rhoAvg = 0.5*(rhoLeft + rhoRight);
    const scalar gammaAvg = 0.5*(gammaLeft + gammaRight);
    const scalar speedOfSoundAvg = sqrt(gammaAvg*pAvg/rhoAvg);
    const scalar alphaEffAvg = 0.5*(alphaEffLeft + alphaEffRight);

    spectralRadiiInv = w_*(mag(contrVavg) + speedOfSoundAvg)*magSf;
    spectralRadiiVis = sqr(magSf)*max(4./3.,gammaAvg)*alphaEffAvg/rhoAvg;

    spectralRadiiown
        = w_*(mag(contrVLeft) + speedOfSoundLeft)*magSf
        + magSf*max(4./3.,gammaLeft)*alphaEffLeft/rhoLeft/mag(dLoc);

    spectralRadiinei
        = w_*(mag(contrVRight) + speedOfSoundRight)*magSf
        + magSf*max(4./3.,gammaRight)*alphaEffRight/rhoRight/mag(dLoc);
}


void spectralRadius::update()
{
    // Get face-to-cell addressing: face area point from owner to neighbour
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    const volVectorField& C = mesh_.C();

    const volScalarField alphaEff(turbulenceModel_.alphaEff());

    const volScalarField gamma
    (
        "gamma",    
        thermoPhysicalModel_.Cp()/thermoPhysicalModel_.Cv()
    );

    // Calculate fluxes at internal faces
    forAll(owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector dLoc = C[nei] - C[own];

        // calculate spectral radius
        spectralRadius::evaluateSpectralRadii
        (
            spectralRadiiInv_[faceI],
            spectralRadiiVis_[faceI],
            spectralRadiiown_[faceI],
            spectralRadiinei_[faceI],

            dLoc,
            alphaEff[own],
            alphaEff[nei],

            p_[own],
            p_[nei],
            U_[own],
            U_[nei],
            rho_[own],
            rho_[nei],

            gamma[own],       // left gamma
            gamma[nei],       // right gamma

            Sf[faceI],      // face vector
            magSf[faceI]   // face area
        ); 
    }

    // Update boundary field and values
    forAll(U_.boundaryField(), patchi)
    {
        const labelUList& pOwner = mesh_.boundary()[patchi].faceCells();

        fvsPatchScalarField& pSpectralRadiiInv
            = spectralRadiiInv_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pSpectralRadiiVis
            = spectralRadiiVis_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pSpectralRadiiown
            = spectralRadiiown_.boundaryFieldRef()[patchi];
        fvsPatchScalarField& pSpectralRadiinei
            = spectralRadiinei_.boundaryFieldRef()[patchi];

        const fvPatchScalarField& palphaEff 
            = alphaEff.boundaryField()[patchi];

        const fvPatchScalarField& pp = p_.boundaryField()[patchi];
        const fvPatchVectorField& pU = U_.boundaryField()[patchi];
        const fvPatchScalarField& prho = rho_.boundaryField()[patchi];

        const fvPatchScalarField& pgamma = gamma.boundaryField()[patchi];

        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchi];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchi];

        // special treatment at coupled boundaries
        if (pp.coupled())
        {
            // primitive variables
            const scalarField palphaEffLeft(palphaEff.patchInternalField());
            const scalarField palphaEffRight(palphaEff.patchNeighbourField());

            const scalarField ppLeft(pp.patchInternalField());
            const scalarField ppRight(pp.patchNeighbourField());

            const vectorField pULeft(pU.patchInternalField());
            const vectorField pURight(pU.patchNeighbourField());

            const scalarField prhoLeft(prho.patchInternalField());
            const scalarField prhoRight(prho.patchNeighbourField());

            const scalarField pgammaLeft(pgamma.patchInternalField());
            const scalarField pgammaRight(pgamma.patchNeighbourField());

            vectorField pd(mesh_.Cf().boundaryField()[patchi].patch().delta());

            forAll(pOwner, faceI)
            {
                // Calculate spectral radius at coupled boundary faces
                spectralRadius::evaluateSpectralRadii
                (
                    pSpectralRadiiInv[faceI],
                    pSpectralRadiiVis[faceI],
                    pSpectralRadiiown[faceI],
                    pSpectralRadiinei[faceI],

                    pd[faceI],
                    palphaEffLeft[faceI],
                    palphaEffRight[faceI],

                    ppLeft[faceI],  // face U
                    ppRight[faceI], // face U
                    pULeft[faceI],  // face U
                    pURight[faceI], // face U
                    prhoLeft[faceI],          // face rho
                    prhoRight[faceI],         // face rho

                    pgammaLeft[faceI],  // face gamma
                    pgammaRight[faceI], // face gamma

                    pSf[faceI],         // face vector
                    pMagSf[faceI]      // face area
                );  
            }
        }
        else
        {
            const vectorField pCf =  pp.patch().Cf();

            forAll(pp, faceI)
            {
                label faceCellI = pp.patch().faceCells()[faceI];
                vector pCLeft = C[faceCellI];

                vector dLoc = (pCf[faceI] - pCLeft);
                // Check if we have to use normal vector at the boundaries ???
                // by Gill

                // Calculate spectral radius at boundary faces
                spectralRadius::evaluateSpectralRadii
                (
                    pSpectralRadiiInv[faceI],
                    pSpectralRadiiVis[faceI],
                    pSpectralRadiiown[faceI],
                    pSpectralRadiinei[faceI],

                    dLoc,
                    alphaEff[faceCellI],
                    alphaEff[faceCellI],

                    p_[faceI],        // face p
                    p_[faceI],        // face p
                    U_[faceI],        // face U
                    U_[faceI],        // face U
                    rho_[faceI],      // face rho
                    rho_[faceI],      // face rho

                    gamma[faceI],    // face gamma
                    gamma[faceI],    // face gamma

                    pSf[faceI],       // face vector
                    pMagSf[faceI]    // face area
                ); 
            }
        }
    }
}


const tmp<scalarField> spectralRadius::LambdaC() const
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    auto tLambdaC = tmp<scalarField>::New(mesh_.V().size(), 0.0);

    scalarField& LambdaC = tLambdaC.ref();

    forAll(owner, facej)
    {
        LambdaC[owner[facej]] += spectralRadiiInv_[facej];
        LambdaC[neighbour[facej]] += spectralRadiiInv_[facej];
    }

    forAll(mesh_.boundary(), patchi)
    {
        const labelUList& pFaceCells 
            = mesh_.boundary()[patchi].faceCells();

        const fvsPatchScalarField& pSpectralRadiiInv
            = spectralRadiiInv_.boundaryField()[patchi];

        forAll(mesh_.boundary()[patchi], facej)
        {
            LambdaC[pFaceCells[facej]] += pSpectralRadiiInv[facej];
        }
    }

    LambdaC = LambdaC/w_;

    return tLambdaC;
}


const tmp<scalarField> spectralRadius::LambdaV() const
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    auto tLambdaV = tmp<scalarField>::New(mesh_.V().size(), 0.0);

    scalarField& LambdaV= tLambdaV.ref();

    forAll(owner, facej)
    {
        LambdaV[owner[facej]] += spectralRadiiVis_[facej];
        LambdaV[neighbour[facej]] += spectralRadiiVis_[facej];
    }

    forAll(mesh_.boundary(), patchi)
    {
        const labelUList& pFaceCells 
            = mesh_.boundary()[patchi].faceCells();

        const fvsPatchScalarField& pSpectralRadiiVis
            = spectralRadiiVis_.boundaryField()[patchi];

        forAll(mesh_.boundary()[patchi], facej)
        {
            LambdaV[pFaceCells[facej]] += pSpectralRadiiVis[facej];
        }
    }

    LambdaV /= mesh_.V();

    return tLambdaV;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
