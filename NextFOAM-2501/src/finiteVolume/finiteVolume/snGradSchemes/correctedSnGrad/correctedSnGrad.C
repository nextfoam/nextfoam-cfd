/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "correctedSnGrad.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "linear.H"
#include "fvcGrad.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fv::correctedSnGrad<Type>::fullGradCorrection
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    if 
    (
        mesh.objectRegistry::foundObject<surfaceScalarField>
        (
            "snGradCorrLimiter"
        )
    )
    {
        typedef typename outerProduct<vector, Type>::type GradType;
        typedef GeometricField<GradType, fvPatchField, volMesh> volGradField;
        typedef GeometricField<GradType, fvsPatchField, surfaceMesh>
            surfaceGradField;

        tmp<volGradField> tgradVf
        (
            new volGradField
            (
                IOobject
                (
                    "tgrad(" + vf.name() + ')',
                    mesh.time().timeName(),
                    mesh
                ),
                mesh,
                dimensioned<GradType>
                (
                    "zero",
                    vf.dimensions()/dimLength,
                    pTraits<GradType>::zero
                )
            )
        );        

        volGradField& gradVf = tgradVf.ref();

        /*if (mesh.foundObject<volGradField>("grad(" + vf.name() + ')'))
        {
            // Use existing gradient field
            gradVf = mesh.lookupObject<volGradField>
            (
                "grad(" + vf.name() + ')'
            );
        }
        else*/
        {
            // Calculate new gradient field
            gradVf = Foam::fv::gradScheme<Type>::New
            (
                mesh,
                IStringStream("VKLimited leastSquares 0.5")()
            )().grad(vf, "grad(" + vf.name() + ')')();
        }

        tmp<surfaceGradField> tgradVsf
        (
            surfaceInterpolationScheme<GradType>::New
            (
                mesh,
                IStringStream("linear")()
            )().interpolate(tgradVf())
        );

        const surfaceScalarField& snGradCorrLimiter
        (
            mesh.objectRegistry::lookupObject<surfaceScalarField>
            (
                "snGradCorrLimiter"
            )
        );

        // construct GeometricField<Type, fvsPatchField, surfaceMesh>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tssf =
            snGradCorrLimiter*mesh.nonOrthCorrectionVectors() & tgradVsf();

        tssf.ref().rename("snGradCorr(" + vf.name() + ')');

        return tssf;
    } // by Gill
    else
    {
        // construct GeometricField<Type, fvsPatchField, surfaceMesh>
        tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tssf =
            linear<typename outerProduct<vector, Type>::type>
            (
                mesh
            ).dotInterpolate
            (
                mesh.nonOrthCorrectionVectors(),
                gradScheme<Type>::New
                (
                    mesh,
                    mesh.gradScheme("grad(" + vf.name() + ')')
                )().grad(vf, "grad(" + vf.name() + ')')
            );
        tssf.ref().rename("snGradCorr(" + vf.name() + ')');

        return tssf;
    }
}


template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh>>
Foam::fv::correctedSnGrad<Type>::correction
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
    const fvMesh& mesh = this->mesh();

    // construct GeometricField<Type, fvsPatchField, surfaceMesh>
    tmp<GeometricField<Type, fvsPatchField, surfaceMesh>> tssf
    (
        new GeometricField<Type, fvsPatchField, surfaceMesh>
        (
            IOobject
            (
                "snGradCorr("+vf.name()+')',
                vf.instance(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            vf.dimensions()*mesh.nonOrthDeltaCoeffs().dimensions()
        )
    );
    GeometricField<Type, fvsPatchField, surfaceMesh>& ssf = tssf.ref();
    ssf.setOriented();

    for (direction cmpt = 0; cmpt < pTraits<Type>::nComponents; ++cmpt)
    {
        ssf.replace
        (
            cmpt,
            correctedSnGrad<typename pTraits<Type>::cmptType>(mesh)
           .fullGradCorrection(vf.component(cmpt))
        );
    }

    return tssf;
}


// ************************************************************************* //
