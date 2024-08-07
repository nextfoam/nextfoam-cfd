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

#include "fvcMeshPhi.H"
#include "fvMesh.H"
#include "ddtScheme.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField> Foam::fvc::meshPhi
(
    const volVectorField& vf
)
{
    return fv::ddtScheme<vector>::New
    (
        vf.mesh(),
        vf.mesh().ddtScheme("ddt(" + vf.name() + ')')
    ).ref().meshPhi(vf);
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvc::meshPhi
(
    const dimensionedScalar& rho,
    const volVectorField& vf
)
{
    return fv::ddtScheme<vector>::New
    (
        vf.mesh(),
        vf.mesh().ddtScheme("ddt(" + rho.name() + ',' + vf.name() + ')')
    ).ref().meshPhi(vf);
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvc::meshPhi
(
    const volScalarField& rho,
    const volVectorField& vf
)
{
    return fv::ddtScheme<vector>::New
    (
        vf.mesh(),
        vf.mesh().ddtScheme("ddt(" + rho.name() + ',' + vf.name() + ')')
    ).ref().meshPhi(vf);
}


void Foam::fvc::makeRelative
(
    surfaceScalarField& phi,
    const volVectorField& U
)
{
    if (phi.mesh().moving())
    {
        phi -= fvc::meshPhi(U);
    }
}

void Foam::fvc::makeRelative
(
    surfaceScalarField& phi,
    const dimensionedScalar& rho,
    const volVectorField& U
)
{
    if (phi.mesh().moving())
    {
        phi -= rho*fvc::meshPhi(rho, U);
    }
}

void Foam::fvc::makeRelative
(
    surfaceScalarField& phi,
    const volScalarField& rho,
    const volVectorField& U
)
{
    if (phi.mesh().moving())
    {
        if (phi.mesh().foundObject<surfaceScalarField>("faceRho"))
        {
            const surfaceScalarField& rhof
            (
                phi.mesh().lookupObject<surfaceScalarField>("faceRho")
            );

            phi -= rhof*fvc::meshPhi(rho, U);
        } // by Gill
        else
        {
            phi -= fvc::interpolate(rho)*fvc::meshPhi(rho, U);
        }
    }
}


void Foam::fvc::makeAbsolute
(
    surfaceScalarField& phi,
    const volVectorField& U
)
{
    if (phi.mesh().moving())
    {
        phi += fvc::meshPhi(U);
    }
}

void Foam::fvc::makeAbsolute
(
    surfaceScalarField& phi,
    const dimensionedScalar& rho,
    const volVectorField& U
)
{
    if (phi.mesh().moving())
    {
        phi += rho*fvc::meshPhi(rho, U);
    }
}

void Foam::fvc::makeAbsolute
(
    surfaceScalarField& phi,
    const volScalarField& rho,
    const volVectorField& U
)
{
    if (phi.mesh().moving())
    {
        if (phi.mesh().foundObject<surfaceScalarField>("faceRho"))
        {
            const surfaceScalarField& rhof
            (
                phi.mesh().lookupObject<surfaceScalarField>("faceRho")
            );

            phi += rhof*fvc::meshPhi(rho, U);
        } // by Gill
        else
        {
            phi += fvc::interpolate(rho)*fvc::meshPhi(rho, U);
        }
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvc::relative
(
    const tmp<surfaceScalarField>& tphi,
    const volVectorField& U
)
{
    if (tphi().mesh().moving())
    {
        return tphi - fvc::meshPhi(U);
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvc::relative
(
    const tmp<surfaceScalarField>& tphi,
    const volScalarField& rho,
    const volVectorField& U
)
{
    if (tphi().mesh().moving())
    {
        if (tphi().mesh().foundObject<surfaceScalarField>("faceRho"))
        {
            const surfaceScalarField& rhof
            (
                tphi().mesh().lookupObject<surfaceScalarField>("faceRho")
            );

            return tphi - rhof*fvc::meshPhi(rho, U);
        } // by Gill

        return tphi - fvc::interpolate(rho)*fvc::meshPhi(rho, U);
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvc::absolute
(
    const tmp<surfaceScalarField>& tphi,
    const volVectorField& U
)
{
    if (tphi().mesh().moving())
    {
        return tphi + fvc::meshPhi(U);
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


Foam::tmp<Foam::surfaceScalarField> Foam::fvc::absolute
(
    const tmp<surfaceScalarField>& tphi,
    const volScalarField& rho,
    const volVectorField& U
)
{
    if (tphi().mesh().moving())
    {
        if (tphi().mesh().foundObject<surfaceScalarField>("faceRho"))
        {
            const surfaceScalarField& rhof
            (
                tphi().mesh().lookupObject<surfaceScalarField>("faceRho")
            );

            return tphi + rhof*fvc::meshPhi(rho, U);
        } // by Gill

        return tphi + fvc::interpolate(rho)*fvc::meshPhi(rho, U);
    }
    else
    {
        return tmp<surfaceScalarField>(tphi, true);
    }
}


void Foam::fvc::correctUf
(
    autoPtr<surfaceVectorField>& Uf,
    const volVectorField& U,
    const surfaceScalarField& phi
)
{
    const fvMesh& mesh = U.mesh();

    if (mesh.dynamic())
    {
        Uf() = fvc::interpolate(U);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        Uf() += n*(phi/mesh.magSf() - (n & Uf()));
    }
}


void Foam::fvc::correctRhoUf
(
    autoPtr<surfaceVectorField>& rhoUf,
    const volScalarField& rho,
    const volVectorField& U,
    const surfaceScalarField& phi
)
{
    const fvMesh& mesh = U.mesh();

    if (mesh.dynamic())
    {
        rhoUf() = fvc::interpolate(rho*U);
        surfaceVectorField n(mesh.Sf()/mesh.magSf());
        rhoUf() += n*(fvc::absolute(phi, rho, U)/mesh.magSf() - (n & rhoUf()));
    }
}


// ************************************************************************* //
