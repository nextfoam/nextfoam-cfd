/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

#include "turbulenceModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "wallFvPatch.H"
#include "fvc.H" // by Gill

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(turbulenceModel, 0);
}

const Foam::word Foam::turbulenceModel::propertiesName("turbulenceProperties");


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulenceModel::turbulenceModel
(
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const word& propertiesName
)
:
    IOdictionary
    (
        IOobject
        (
            IOobject::groupName(propertiesName, alphaRhoPhi.group()),
            U.time().constant(),
            U.db(),
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    ),

    runTime_(U.time()),
    mesh_(U.mesh()),

    U_(U),
    alphaRhoPhi_(alphaRhoPhi),
    phi_(phi),
    y_(mesh_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::turbulenceModel::phi() const
{
    return phi_;
}


Foam::tmp<Foam::volScalarField::Internal> Foam::turbulenceModel::divU() const
{
    tmp<volScalarField::Internal> tDivU
    (
        new volScalarField::Internal
        (
            turbulenceModel::phi()().dimensions() == dimMass/dimTime
          ? fvc::div
            (
                fvc::absolute
                (
                    turbulenceModel::phi(), 
                    mesh().lookupObject<volScalarField>("rho"), 
                    U()
                )
            )().v()
          : fvc::div(fvc::absolute(turbulenceModel::phi(), U()))().v()
        )
    ); 

    volScalarField::Internal& divU = tDivU.ref();

    if (turbulenceModel::phi()().dimensions() == dimMass/dimTime)
    {
        const volScalarField& rho
        (   
            mesh().lookupObject<volScalarField>("rho")
        );

        divU -= (U() & fvc::grad(rho))().v();
        
        divU /= rho.v();
    }

    return tDivU;
} // by Gill


bool Foam::turbulenceModel::read()
{
    return regIOobject::read();
}


void Foam::turbulenceModel::validate()
{}


void Foam::turbulenceModel::correct()
{
    if (mesh_.changing())
    {
        y_.correct();
    }
}


// ************************************************************************* //
