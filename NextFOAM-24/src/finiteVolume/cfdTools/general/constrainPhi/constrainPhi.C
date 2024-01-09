/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2016 OpenFOAM Foundation
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

#include "constrainPhi.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "fixedFluxExtrapolatedPressureFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::constrainPhi
(
    const tmp<surfaceScalarField>& tPhi,
    const volVectorField& U,
    const volScalarField& p
)
{
    tmp<surfaceScalarField> tPhiNew;

    if (tPhi.isTmp())
    {
        tPhiNew = tPhi;
        tPhiNew.ref().rename("Phi");
    }
    else
    {
        tPhiNew = new surfaceScalarField("Phi", tPhi);
    }

    surfaceScalarField& Phi = tPhiNew.ref();
    surfaceScalarField::Boundary& Phibf = Phi.boundaryFieldRef();

    forAll(U.boundaryField(), patchi)
    {
        if
        (
           !U.boundaryField()[patchi].assignable()
        && !isA<fixedFluxExtrapolatedPressureFvPatchScalarField>
            (
                p.boundaryField()[patchi]
            )
        )
        {
            Phibf[patchi] = 0.0;
        }
    }

    return tPhiNew;
}


Foam::tmp<Foam::surfaceScalarField> Foam::constrainPhig
(
    const tmp<surfaceScalarField>& tPhig,
    const volVectorField& U
)
{
    tmp<surfaceScalarField> tPhigNew;

    if (tPhig.isTmp())
    {
        tPhigNew = tPhig;
        tPhigNew.ref().rename("Phig");
    }
    else
    {
        tPhigNew = new surfaceScalarField("Phig", tPhig);
    }

    surfaceScalarField& Phig = tPhigNew.ref();
    surfaceScalarField::Boundary& Phigbf = Phig.boundaryFieldRef();

    forAll(U.boundaryField(), patchi)
    {
        if
        (
           !U.boundaryField()[patchi].assignable()
        && (U.boundaryField()[patchi].patch().patch().type() == "patch")
        )
        {
            Phigbf[patchi] = 0.0;
        }
    }

    return tPhigNew;
}


// ************************************************************************* //
