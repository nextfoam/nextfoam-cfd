/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2019-2023 OpenCFD Ltd.
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

#include "emptyFaPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void Foam::fa::jouleHeatingSource::initialiseSigma
(
    const dictionary& dict,
    autoPtr<Function1<Type>>& sigmaVsTPtr
)
{
    typedef GeometricField<Type, faPatchField, areaMesh> AreaFieldType;

    auto& obr = regionMesh().thisDb();

    if (dict.found("sigma"))
    {
        // Sigma to be defined using a Function1 type
        sigmaVsTPtr = Function1<Type>::New("sigma", dict, &mesh_);

        auto tsigma = tmp<AreaFieldType>::New
        (
            IOobject
            (
                typeName + ":sigma_" + regionName_,
                obr.time().timeName(),
                obr,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            regionMesh(),
            dimensioned<Type>(sqr(dimCurrent)/dimPower/dimLength, Zero)
        );

        regIOobject::store(tsigma.ptr());

        Info<< "    Conductivity 'sigma' read from dictionary as f(T)"
            << nl << endl;
    }
    else
    {
        // Sigma to be defined by user input
        auto tsigma = tmp<AreaFieldType>::New
        (
            IOobject
            (
                typeName + ":sigma_" + regionName_,
                obr.time().timeName(),
                obr,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            regionMesh()
        );

        regIOobject::store(tsigma.ptr());

        Info<< "    Conductivity 'sigma' read from file" << nl << endl;
    }
}


template<class Type>
const Foam::GeometricField<Type, Foam::faPatchField, Foam::areaMesh>&
Foam::fa::jouleHeatingSource::updateSigma
(
    const autoPtr<Function1<Type>>& sigmaVsTPtr
) const
{
    typedef GeometricField<Type, faPatchField, areaMesh> AreaFieldType;

    const auto& obr = regionMesh().thisDb();

    auto& sigma =
        obr.lookupObjectRef<AreaFieldType>(typeName + ":sigma_" + regionName_);

    if (!sigmaVsTPtr)
    {
        // Electrical conductivity field, sigma, was specified by the user
        return sigma;
    }

    const auto& T = obr.lookupObject<areaScalarField>(TName_);

    // Internal field
    forAll(sigma, i)
    {
        sigma[i] = sigmaVsTPtr->value(T[i]);
    }


    // Boundary field
    typename AreaFieldType::Boundary& bf = sigma.boundaryFieldRef();
    forAll(bf, patchi)
    {
        faPatchField<Type>& pf = bf[patchi];
        if (!isA<emptyFaPatch>(bf[patchi]))
        {
            const scalarField& Tbf = T.boundaryField()[patchi];
            forAll(pf, facei)
            {
                pf[facei] = sigmaVsTPtr->value(Tbf[facei]);
            }
        }
    }

    // Update processor patches
    sigma.correctBoundaryConditions();

    return sigma;
}


// ************************************************************************* //
