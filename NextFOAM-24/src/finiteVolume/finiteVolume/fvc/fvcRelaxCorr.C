/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "fvcRelaxCorr.H"
#include "fvMesh.H"
#include "momentum.H"
#include "fixedFluxExtrapolatedPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::fvc::relaxCorr
(
    const surfaceScalarField& phi,
    const volVectorField& U,
    const volScalarField& p,
    const word& name
)
{
    const fvMesh& mesh = U.mesh();

    tmp<surfaceScalarField> tprevUPhi
    (   
        momentumInterpolate(U.prevIter()) & mesh.Sf()
    );

    tmp<surfaceScalarField> tPhiCorr
    (
        new surfaceScalarField
        (
            IOobject
            (
                "relaxCorr(" + U.name() + ',' + phi.name() + ')',
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensioned<scalar>
            (   
                "zero", 
                phi.dimensions(), 
                pTraits<scalar>::zero
            )
        )
    );

    scalarField& phiCorr
    (   
        tPhiCorr.ref().primitiveFieldRef()
    );

    phiCorr = phi.prevIter().primitiveField() - tprevUPhi().primitiveField();

    surfaceScalarField::Boundary& bPhiCorr = tPhiCorr.ref().boundaryFieldRef();

    forAll(U.boundaryField(), patchi)
    {
        if  
        (   
            U.boundaryField()[patchi].assignable()
        ||  isA<fixedFluxExtrapolatedPressureFvPatchScalarField>
            (   
                p.boundaryField()[patchi]
            )
        )
        {
            bPhiCorr[patchi] =
                phi.prevIter().boundaryField()[patchi]
              - tprevUPhi().boundaryField()[patchi];
        }
    }

    scalar UUrf(1.0);

    if (mesh.relaxEquation(name))
    {
        UUrf = mesh.equationRelaxationFactor(name);
    }

    return (1.0 - UUrf)*tPhiCorr;
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvc::relaxCorr
(
    const surfaceScalarField& rhof,
    const surfaceScalarField& phi,
    const volVectorField& U,
    const volScalarField& p,
    const word& name
)
{
    const fvMesh& mesh = U.mesh();

    tmp<surfaceScalarField> tprevUPhi
    (
        rhof*(momentumInterpolate(U.prevIter()) & mesh.Sf())
    );

    tmp<surfaceScalarField> tPhiCorr
    (
        new surfaceScalarField
        (
            IOobject
            (
                "relaxCorr(" + U.name() + ',' + phi.name() + ')',
                mesh.time().timeName(),
                mesh
            ),
            mesh,
            dimensioned<scalar>
            (
                "zero",
                phi.dimensions(),
                pTraits<scalar>::zero
            )
        )
    );

    scalarField& phiCorr
    (
        tPhiCorr.ref().primitiveFieldRef()
    );

    phiCorr = phi.prevIter().primitiveField() - tprevUPhi().primitiveField();

    surfaceScalarField::Boundary& bPhiCorr = tPhiCorr.ref().boundaryFieldRef();

    forAll(U.boundaryField(), patchi)
    {
        if
        (
            U.boundaryField()[patchi].assignable()
        ||  isA<fixedFluxExtrapolatedPressureFvPatchScalarField>
            (
                p.boundaryField()[patchi]
            )
        )
        {
            bPhiCorr[patchi] =
                phi.prevIter().boundaryField()[patchi]
              - tprevUPhi().boundaryField()[patchi];
        }
    }

    scalar UUrf(1.0);

    if (mesh.relaxEquation(name))
    {
        UUrf = mesh.equationRelaxationFactor(name);
    }

    return (1.0 - UUrf)*tPhiCorr;
}


// ************************************************************************* //
