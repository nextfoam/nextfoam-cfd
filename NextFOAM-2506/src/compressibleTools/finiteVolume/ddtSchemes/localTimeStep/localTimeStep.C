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

#include "localTimeStep.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace compressible
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

localTimeStep::localTimeStep(const volVectorField& U)
:
    mesh_(U.mesh()),
    U_(U),
    CoDeltaT_
    (
        IOobject
        (
            "CoDeltaT",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("CoDeltaT", dimTime, 0.1),
        fvPatchFieldBase::zeroGradientType()
    ),
    trDeltaT_
    (
        new volScalarField
        (
            IOobject
            (
                "rDeltaT",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("small", dimless/dimTime, SMALL),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    ),
    deltaS_
    (
        IOobject
        (
            "deltaS",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("deltaS", dimArea, vector::one),
        fvPatchFieldBase::zeroGradientType()
    ),
    cellVolume_
    (
        IOobject
        (
            "cellVolume",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("cellVolume", dimVolume, 0.1),
        fvPatchFieldBase::zeroGradientType()
    )
{
    updateDeltaS();
};


void localTimeStep::updateDeltaS()
{
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    // new formulation for deltaS according to Blazek

    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();

    // Reset values
    deltaS_ = dimensionedVector("deltaS", dimArea, vector::zero);
    cellVolume_.primitiveFieldRef() = mesh_.V();
    cellVolume_.correctBoundaryConditions();

    forAll(owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        deltaS_[own] += 0.5*cmptMag(Sf[faceI]);
        deltaS_[nei] += 0.5*cmptMag(Sf[faceI]);
    }

    forAll(deltaS_.boundaryField(), patchI)
    {
        const fvPatchVectorField& pp = deltaS_.boundaryField()[patchI];

        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchI];

        const labelUList& faceCells =  pp.patch().faceCells();

        forAll(pp, faceI)
        {
            label own = faceCells[faceI];

            deltaS_[own] += 0.5*cmptMag(pSf[faceI]);
        }
    }

    deltaS_.correctBoundaryConditions();
};


void localTimeStep::update
(
    scalar maxCo, 
    const tmp<scalarField> Lambda, 
    Switch adjustLocalTimeStep
)
{
    if (mesh_.moving())
    {
        updateDeltaS();
    }
    
    CoDeltaT_.primitiveFieldRef() = maxCo*mesh_.V()/Lambda;

    if (adjustLocalTimeStep)
    {
        CoDeltaT_.primitiveFieldRef() = min(CoDeltaT_.primitiveField());
    }

    CoDeltaT_.correctBoundaryConditions();

    volScalarField& rDeltaT = trDeltaT_.ref();

    scalar rDeltaTSmoothingCoeff
    (
        mesh_.time().controlDict().lookupOrDefault<scalar>
        (
            "rDeltaTSmoothingCoeff",
            1.0
        )
    );

    scalar rDeltaTDampingCoeff
    (
        mesh_.time().controlDict().lookupOrDefault<scalar>
        (
            "rDeltaTDampingCoeff",
            1.0
        )
    );

    volScalarField rDeltaT0("rDeltaT0", rDeltaT);

    rDeltaT.primitiveFieldRef() =
        1.0/CoDeltaT_.primitiveField();

    // Update tho boundary values of the reciprocal time-step
    rDeltaT.correctBoundaryConditions();

    Info<< "Flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    if (rDeltaTSmoothingCoeff < 1.0)
    {
        fvc::smooth(rDeltaT, rDeltaTSmoothingCoeff);
    }

    Info<< "Smoothed flow time scale min/max = "
        << gMin(1/rDeltaT.primitiveField())
        << ", " << gMax(1/rDeltaT.primitiveField()) << endl;

    // Limit rate of change of time scale
    // - reduce as much as required
    // - only increase at a fraction of old time scale
    if
    (
        rDeltaTDampingCoeff < 1.0
     && mesh_.time().timeIndex() > mesh_.time().startTimeIndex() + 1
    )
    {
        rDeltaT =
            rDeltaT0
           *max(rDeltaT/rDeltaT0, scalar(1) - rDeltaTDampingCoeff);

        Info<< "Damped flow time scale min/max = "
            << gMin(1/rDeltaT.primitiveField())
            << ", " << gMax(1/rDeltaT.primitiveField()) << endl;
    }

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
