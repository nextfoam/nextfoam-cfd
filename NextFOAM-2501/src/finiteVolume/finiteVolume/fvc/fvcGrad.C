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

#include "fvcGrad.H"
#include "fvcSurfaceIntegrate.H"
#include "fvMesh.H"
#include "gaussGrad.H"
#include "fvcReconstruct.H"
#include "fvcSnGrad.H"
#include "linearPressure.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fvc
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >
>
grad
(
    const GeometricField<Type, fvsPatchField, surfaceMesh>& ssf
)
{
    return fv::gaussGrad<Type>::gradf(ssf, "grad(" + ssf.name() + ')');
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, fvPatchField, volMesh
    >
>
grad
(
    const tmp<GeometricField<Type, fvsPatchField, surfaceMesh>>& tssf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    tmp<GeometricField<GradType, fvPatchField, volMesh>> Grad
    (
        fvc::grad(tssf())
    );
    tssf.clear();
    return Grad;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, fvPatchField, volMesh
    >
>
grad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    return fv::gradScheme<Type>::New
    (
        vf.mesh(),
        vf.mesh().gradScheme(name)
    )().grad(vf, name);
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, fvPatchField, volMesh
    >
>
grad
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf,
    const word& name
)
{
    tmp
    <
        GeometricField
        <
            typename outerProduct<vector, Type>::type, fvPatchField, volMesh
        >
    > tGrad
    (
        fvc::grad(tvf(), name)
    );
    tvf.clear();
    return tGrad;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, fvPatchField, volMesh
    >
>
grad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    return fvc::grad(vf, "grad(" + vf.name() + ')');
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, fvPatchField, volMesh
    >
>
grad
(
    const tmp<GeometricField<Type, fvPatchField, volMesh>>& tvf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    tmp<GeometricField<GradType, fvPatchField, volMesh>> Grad
    (
        fvc::grad(tvf())
    );
    tvf.clear();
    return Grad;
}


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, fvPatchField, volMesh
    >
>
calcGaussGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf,
    const word& name
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh>> tgGrad
    (
        div(interpolate(vf) * mesh.Sf())
    );

    GeometricField<GradType, fvPatchField, volMesh>& gGrad = tgGrad.ref();

    gGrad.rename(name);

    gGrad.correctBoundaryConditions();

    typename GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >::Boundary& gGradbf = gGrad.boundaryFieldRef();

    forAll(vf.boundaryField(), patchi)
    {
        if (!vf.boundaryField()[patchi].coupled())
        {
            const vectorField n
            (
                mesh.Sf().boundaryField()[patchi]
              / mesh.magSf().boundaryField()[patchi]
            );

            gGradbf[patchi] += n *
            (
                vf.boundaryField()[patchi].snGrad()
              - (n & gGradbf[patchi])
            );
        }
    }

    return tgGrad;
} // by Gill


template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, fvPatchField, volMesh
    >
>
gaussGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;
    typedef GeometricField<GradType, fvPatchField, volMesh> GradFieldType;

    const word name("grad(" + vf.name() +")");

    GradFieldType* pgGrad =
        vf.mesh().objectRegistry::template getObjectPtr<GradFieldType>(name);

    if (!vf.mesh().cache(name) || vf.mesh().changing())
    {
        // Delete any old occurrences to avoid double registration
        if (pgGrad && pgGrad->ownedByRegistry())
        {
            solution::cachePrintMessage("Deleting", name, vf);
            delete pgGrad;
        }

        solution::cachePrintMessage("Calculating", name, vf);
        return calcGaussGrad(vf, name);
    }


    if (!pgGrad)
    {
        solution::cachePrintMessage("Calculating and caching", name, vf);

        pgGrad = calcGaussGrad(vf, name).ptr();
        regIOobject::store(pgGrad);
    }
    else
    {
        if (pgGrad->upToDate(vf))
        {
            solution::cachePrintMessage("Reusing", name, vf);
        }
        else
        {
            solution::cachePrintMessage("Updating", name, vf);
            delete pgGrad;

            pgGrad = calcGaussGrad(vf, name).ptr();
            regIOobject::store(pgGrad);
        }
    }

    return *pgGrad;
} // by Gill


//- This kind of pressure gradient evaluation method is not neeed if we find 
//- suitable pressure interpolation method for gauss grad scheme when gravity 
//- is introduced
template<class Type>
tmp
<
    GeometricField
    <
        typename outerProduct<vector,Type>::type, fvPatchField, volMesh
    >
>
reconSnGrad
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    typedef typename outerProduct<vector, Type>::type GradType;

    const fvMesh& mesh = vf.mesh();

    tmp<GeometricField<GradType, fvPatchField, volMesh>> tsGrad
    (
        fvc::reconstruct(fvc::snGrad(vf) * mesh.magSf())
    );

    GeometricField<GradType, fvPatchField, volMesh>& sGrad = tsGrad.ref();

    // Adjust pressure gradient at porous cells
    if (mesh.foundObject<fv::optionList>("fvOptions"))
    {
        const fv::optionList& fvOptions
        (
            mesh.lookupObject<fv::optionList>("fvOptions")
        );

        if (fvOptions.hasActivePorousZone())
        {
            tmp<surfaceScalarField> tLinearFacePressure
            (
                pressureInterpolationScheme::interpolate
                (
                    vf,
                    linearPressure(mesh).weights(vf)
                )
            );

            tLinearFacePressure.ref() += linearPressure(mesh).correction(vf);

            // - calculate gaussGrad only at porous media zone interface cells
            const cellList& cells = mesh.cells();
            const labelUList& owner = mesh.owner();
            const surfaceVectorField& Sf = mesh.Sf();
            const surfaceVectorField::Boundary& bSf = Sf.boundaryField();
            const surfaceScalarField& Lp = tLinearFacePressure();
            const surfaceScalarField::Boundary& bLp = Lp.boundaryField();

            const labelList& porousCells = fvOptions.activePorousCells();

            forAll(porousCells, i)
            {
                vector porousCellGrad(Zero);

                const label& celli(porousCells[i]);

                const cell& cell = cells[celli];

                forAll(cell, j)
                {
                    const label& facei = cell[j];   
                    if (mesh.isInternalFace(facei))
                    {
                        if (celli == owner[facei])
                        {
                            porousCellGrad += (Lp[facei] * Sf[facei]);
                        }
                        else
                        {
                            porousCellGrad -= (Lp[facei] * Sf[facei]);
                        }
                    }
                    else
                    {
                        const label& patchI 
                            = mesh.boundaryMesh().whichPatch(facei);

                        const polyPatch& p = mesh.boundary()[patchI].patch();

                        if (!isA<emptyPolyPatch>(p))
                        {
                            const label& pFaceI
                                = mesh.boundaryMesh()[patchI].whichFace(facei);

                            porousCellGrad
                                += (bLp[patchI][pFaceI] * bSf[patchI][pFaceI]);
                        }
                    }
                }

                porousCellGrad /= mesh.Vsc()()[celli];

                sGrad[celli] = porousCellGrad;
            }
        }
    }

    sGrad.correctBoundaryConditions();

    typename GeometricField
    <
        typename outerProduct<vector, Type>::type, fvPatchField, volMesh
    >::Boundary& sGradbf = sGrad.boundaryFieldRef();

    forAll(vf.boundaryField(), patchi)
    {
        if (!vf.boundaryField()[patchi].coupled())
        {
            const vectorField n
            (
                mesh.Sf().boundaryField()[patchi]
              / mesh.magSf().boundaryField()[patchi]
            );

            sGradbf[patchi] += n *
            (
                vf.boundaryField()[patchi].snGrad()
              - (n & sGradbf[patchi])
            );
        }
    }

    return tsGrad;
} // by Gill


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fvc

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
