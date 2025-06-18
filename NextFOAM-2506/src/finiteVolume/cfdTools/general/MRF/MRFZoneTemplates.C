/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "MRFZone.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::MRFZone::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    surfaceScalarField& phi
) const
{
    if (!active_)
    {
        return;
    }

    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const vector Omega = this->Omega();
    const vector Ut = this->Ut();

    const vectorField& Cfi = Cf;
    const vectorField& Sfi = Sf;
    scalarField& phii = phi.primitiveFieldRef();

    // Internal faces
    forAll(internalFaces_, i)
    {
        label facei = internalFaces_[i];

        vector frameU = Ut + (Omega ^ (Cfi[facei] - origin_));

        phii[facei] -= rho[facei]*frameU & Sfi[facei];
    }

    makeRelativeRhoFlux(rho.boundaryField(), phi.boundaryFieldRef());
}


template<class RhoFieldType>
void Foam::MRFZone::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    FieldField<fvsPatchField, scalar>& phi
) const
{
    if (!active_)
    {
        return;
    }

    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const vector Omega = this->Omega();
    const vector Ut = this->Ut();

    // Included patches
    forAll(includedFaces_, patchi)
    {
        forAll(includedFaces_[patchi], i)
        {
            label patchFacei = includedFaces_[patchi][i];

            phi[patchi][patchFacei] = 0.0;
        }
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        forAll(excludedFaces_[patchi], i)
        {
            label patchFacei = excludedFaces_[patchi][i];

            const vector& Cfb = Cf.boundaryField()[patchi][patchFacei];
            const vector& Sfb = Sf.boundaryField()[patchi][patchFacei];

            vector frameU = Ut + (Omega ^ (Cfb - origin_));

            phi[patchi][patchFacei] -= rho[patchi][patchFacei]*frameU & Sfb;
        }
    }
}


template<class RhoFieldType>
void Foam::MRFZone::makeRelativeRhoFlux
(
    const RhoFieldType& rho,
    Field<scalar>& phi,
    const label patchi
) const
{
    if (!active_)
    {
        return;
    }

    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const vector Omega = this->Omega();
    const vector Ut = this->Ut();

    // Included patches
    forAll(includedFaces_[patchi], i)
    {
        label patchFacei = includedFaces_[patchi][i];

        phi[patchFacei] = 0.0;
    }

    // Excluded patches
    forAll(excludedFaces_[patchi], i)
    {
        label patchFacei = excludedFaces_[patchi][i];

        const vector& Cfb = Cf.boundaryField()[patchi][patchFacei];
        const vector& Sfb = Sf.boundaryField()[patchi][patchFacei];

        vector frameU = Ut + (Omega ^ (Cfb - origin_));

        phi[patchFacei] -= rho[patchFacei]*frameU & Sfb;
    }
}


template<class RhoFieldType>
void Foam::MRFZone::makeAbsoluteRhoFlux
(
    const RhoFieldType& rho,
    surfaceScalarField& phi
) const
{
    if (!active_)
    {
        return;
    }

    const surfaceVectorField& Cf = mesh_.Cf();
    const surfaceVectorField& Sf = mesh_.Sf();

    const vector Omega = this->Omega();
    const vector Ut = this->Ut();

    const vectorField& Cfi = Cf;
    const vectorField& Sfi = Sf;
    scalarField& phii = phi.primitiveFieldRef();

    // Internal faces
    forAll(internalFaces_, i)
    {
        label facei = internalFaces_[i];

        vector frameU = Ut + (Omega ^ (Cfi[facei] - origin_));

        phii[facei] += rho[facei]*frameU & Sfi[facei];
    }

    surfaceScalarField::Boundary& phibf = phi.boundaryFieldRef();


    // Included patches
    forAll(includedFaces_, patchi)
    {
        forAll(includedFaces_[patchi], i)
        {
            label patchFacei = includedFaces_[patchi][i];

            const vector& Cfb = Cf.boundaryField()[patchi][patchFacei];
            const vector& Sfb = Sf.boundaryField()[patchi][patchFacei];

            vector frameU = Ut + (Omega ^ (Cfb - origin_));

            phibf[patchi][patchFacei] +=
                rho.boundaryField()[patchi][patchFacei]*frameU & Sfb;
        }
    }

    // Excluded patches
    forAll(excludedFaces_, patchi)
    {
        forAll(excludedFaces_[patchi], i)
        {
            label patchFacei = excludedFaces_[patchi][i];

            const vector& Cfb = Cf.boundaryField()[patchi][patchFacei];
            const vector& Sfb = Sf.boundaryField()[patchi][patchFacei];

            vector frameU = Ut + (Omega ^ (Cfb - origin_));

            phibf[patchi][patchFacei] +=
                rho.boundaryField()[patchi][patchFacei]*frameU & Sfb;
        }
    }
}


template<class Type>
void Foam::MRFZone::zero
(
    GeometricField<Type, fvsPatchField, surfaceMesh>& phi
) const
{
    if (!active_)
    {
        return;
    }

    Field<Type>& phii = phi.primitiveFieldRef();

    forAll(internalFaces_, i)
    {
        phii[internalFaces_[i]] = Zero;
    }

    auto& phibf = phi.boundaryFieldRef();

    forAll(includedFaces_, patchi)
    {
        forAll(includedFaces_[patchi], i)
        {
            phibf[patchi][includedFaces_[patchi][i]] = Zero;
        }
    }

    forAll(excludedFaces_, patchi)
    {
        forAll(excludedFaces_[patchi], i)
        {
            phibf[patchi][excludedFaces_[patchi][i]] = Zero;
        }
    }
}


// ************************************************************************* //
