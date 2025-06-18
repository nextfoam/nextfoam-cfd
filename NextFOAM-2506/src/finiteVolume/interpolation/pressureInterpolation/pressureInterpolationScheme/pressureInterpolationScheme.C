/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2021 OpenCFD Ltd.
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

#include "pressureInterpolationScheme.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "geometricOneField.H"
#include "coupledFvPatchField.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(pressureInterpolationScheme, 0);
    defineRunTimeSelectionTable(pressureInterpolationScheme, Mesh);
}


// * * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * //

Foam::tmp<Foam::pressureInterpolationScheme>
Foam::pressureInterpolationScheme::New
(
    const fvMesh& mesh,
    Istream& schemeData
)
{
    if (schemeData.eof())
    {
        FatalIOErrorInFunction(schemeData)
            << "Discretisation scheme not specified\n\n"
            << "Valid schemes:\n"
            << MeshConstructorTablePtr_->sortedToc()
            << exit(FatalIOError);
    }

    const word schemeName(schemeData);

    if (pressureInterpolationScheme::debug)
    {
        InfoInFunction << "Discretisation scheme = " << schemeName << endl;
    }

    auto* ctorPtr = MeshConstructorTable(schemeName);

    if (!ctorPtr)
    {
        FatalIOErrorInLookup
        (
            schemeData,
            "discretisation",
            schemeName,
            *MeshConstructorTablePtr_
        ) << exit(FatalIOError);
    }

    return ctorPtr(mesh, schemeData);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::surfaceScalarField>
Foam::pressureInterpolationScheme::interpolate
(
    const volScalarField& vf,
    const tmp<surfaceScalarField>& tlambdas
)
{
    if (pressureInterpolationScheme::debug)
    {
        InfoInFunction
            << "Interpolating "
            << vf.name()
            << " from cells to faces without explicit correction"
            << endl;
    }

    const surfaceScalarField& lambdas = tlambdas();

    const scalarField& vfi = vf;
    const scalarField& lambda = lambdas;

    const fvMesh& mesh = vf.mesh();
    const labelUList& P = mesh.owner();
    const labelUList& N = mesh.neighbour();

    tmp<GeometricField<scalar, fvsPatchField, surfaceMesh>> tsf
    (
        new GeometricField<scalar, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "interpolate("+vf.name()+')',
                vf.instance(),
                vf.db()
            ),
            mesh,
            vf.dimensions()
        )
    );
    GeometricField<scalar, fvsPatchField, surfaceMesh>& sf = tsf.ref();

    Field<scalar>& sfi = sf.primitiveFieldRef();

    for (label fi=0; fi<P.size(); fi++)
    {
        sfi[fi] = lambda[fi]*(vfi[P[fi]] - vfi[N[fi]]) + vfi[N[fi]];
    }

    // Interpolate across coupled patches using given lambdas

    surfaceScalarField::Boundary& sfbf = sf.boundaryFieldRef();

    forAll(lambdas.boundaryField(), pi)
    {
        const fvsPatchScalarField& pLambda = lambdas.boundaryField()[pi];
        fvsPatchField<scalar>& psf = sfbf[pi];

        if (vf.boundaryField()[pi].coupled())
        {
            psf =
                pLambda*vf.boundaryField()[pi].patchInternalField()
              + (1.0 - pLambda)*vf.boundaryField()[pi].patchNeighbourField();
        }
        else
        {
            psf = vf.boundaryField()[pi];
        }
    }

    tlambdas.clear();

    return tsf;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::pressureInterpolationScheme::interpolate(const volScalarField& vf) const
{
    if (pressureInterpolationScheme::debug)
    {
        InfoInFunction
            << "Interpolating "
            << vf.name()
            << " from cells to faces"
            << endl;
    }

    tmp<surfaceScalarField> tsf = interpolate(vf, weights(vf));

    if (corrected())
    {
        tsf.ref() += correction(vf);
    }

    return tsf;
}


Foam::tmp<Foam::surfaceScalarField>
Foam::pressureInterpolationScheme::interpolate
(
    const tmp<volScalarField>& tvf
) const
{
    tmp<surfaceScalarField> tinterpVf = interpolate(tvf());
    tvf.clear();
    return tinterpVf;
}


// ************************************************************************* //
