/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2017,2023 OpenCFD Ltd.
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

#include "backwardDualTimeDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
scalar backwardDualTimeDdtScheme<Type>::deltaT_() const
{
    return mesh().time().deltaTValue();
}


template<class Type>
scalar backwardDualTimeDdtScheme<Type>::deltaT0_() const
{
    return mesh().time().deltaT0Value();
}


template<class Type>
template<class GeoField>
scalar backwardDualTimeDdtScheme<Type>::deltaT0_(const GeoField& vf) const
{
    if (mesh().time().timeIndex() < 2)
    {
        return GREAT;
    }
    else
    {
        return deltaT0_();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
backwardDualTimeDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_();

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                dimensioned<Type>(dt.dimensions()/dimTime, Zero)
            )
        );

        tdtdt.ref().primitiveFieldRef() = rDeltaT.value()*dt.value()*
        (
            coefft - (coefft0*mesh().V0() - coefft00*mesh().V00())/mesh().V()
        );

        return tdtdt;
    }

    return tmp<GeometricField<Type, fvPatchField, volMesh>>::New
    (
        ddtIOobject,
        mesh(),
        dimensioned<Type>(dt.dimensions()/dimTime, Zero),
        fvPatchFieldBase::calculatedType()
    );
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
backwardDualTimeDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    coefft*vf.primitiveField() -
                    (
                        coefft0*vf.oldTime().primitiveField()*mesh().V0()
                      - coefft00*vf.oldTime().oldTime().primitiveField()
                       *mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*
                (
                    coefft*vf.boundaryField() -
                    (
                        coefft0*vf.oldTime().boundaryField()
                      - coefft00*vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );

        // Different operation on boundary v.s. internal so re-evaluate
        // coupled boundaries
        tdtdt.ref().boundaryFieldRef().
            template evaluateCoupled<coupledFvPatch>();

        return tdtdt;
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                  + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
backwardDualTimeDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*rho.value()*
                (
                    coefft*vf.primitiveField() -
                    (
                        coefft0*vf.oldTime().primitiveField()*mesh().V0()
                      - coefft00*vf.oldTime().oldTime().primitiveField()
                       *mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*rho.value()*
                (
                    coefft*vf.boundaryField() -
                    (
                        coefft0*vf.oldTime().boundaryField()
                      - coefft00*vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );

        // Different operation on boundary v.s. internal so re-evaluate
        // coupled boundaries
        tdtdt.ref().boundaryFieldRef().
            template evaluateCoupled<coupledFvPatch>();

        return tdtdt;
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*rho*
                (
                    coefft*vf
                  - coefft0*vf.oldTime()
                 + coefft00*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
backwardDualTimeDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    coefft*rho.primitiveField()*vf.primitiveField() -
                    (
                        coefft0*rho.oldTime().primitiveField()
                       *vf.oldTime().primitiveField()*mesh().V0()
                      - coefft00*rho.oldTime().oldTime().primitiveField()
                       *vf.oldTime().oldTime().primitiveField()*mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*
                (
                    coefft*rho.boundaryField()*vf.boundaryField() -
                    (
                        coefft0*rho.oldTime().boundaryField()
                       *vf.oldTime().boundaryField()
                      - coefft00*rho.oldTime().oldTime().boundaryField()
                       *vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );

        // Different operation on boundary v.s. internal so re-evaluate
        // coupled boundaries
        tdtdt.ref().boundaryFieldRef().
            template evaluateCoupled<coupledFvPatch>();

        return tdtdt;
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*
                (
                    coefft*rho*vf
                  - coefft0*rho.oldTime()*vf.oldTime()
                  + coefft00*rho.oldTime().oldTime()*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
backwardDualTimeDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

    IOobject ddtIOobject
    (
        "ddt("+alpha.name()+','+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    if (mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()
               *alpha.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.value()*
                (
                    coefft
                   *alpha.primitiveField()
                   *rho.primitiveField()
                   *vf.primitiveField() -
                    (
                        coefft0
                       *alpha.oldTime().primitiveField()
                       *rho.oldTime().primitiveField()
                       *vf.oldTime().primitiveField()*mesh().V0()

                      - coefft00
                       *alpha.oldTime().oldTime().primitiveField()
                       *rho.oldTime().oldTime().primitiveField()
                       *vf.oldTime().oldTime().primitiveField()*mesh().V00()
                    )/mesh().V()
                ),
                rDeltaT.value()*
                (
                    coefft
                   *alpha.boundaryField()
                   *rho.boundaryField()
                   *vf.boundaryField() -
                    (
                        coefft0
                       *alpha.oldTime().boundaryField()
                       *rho.oldTime().boundaryField()
                       *vf.oldTime().boundaryField()

                      - coefft00
                       *alpha.oldTime().oldTime().boundaryField()
                       *rho.oldTime().oldTime().boundaryField()
                       *vf.oldTime().oldTime().boundaryField()
                    )
                )
            )
        );

        // Different operation on boundary v.s. internal so re-evaluate
        // coupled boundaries
        tdtdt.ref().boundaryFieldRef().
            template evaluateCoupled<coupledFvPatch>();

        return tdtdt;
    }
    else
    {
        return tmp<GeometricField<Type, fvPatchField, volMesh>>
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                rDeltaT*
                (
                    coefft*alpha*rho*vf
                  - coefft0*alpha.oldTime()*rho.oldTime()*vf.oldTime()
                  + coefft00*alpha.oldTime().oldTime()
                   *rho.oldTime().oldTime()*vf.oldTime().oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type>>
backwardDualTimeDdtScheme<Type>::fvmDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            vf.dimensions()*dimVol/dimTime
        )
    );

    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    const volScalarField& rDeltaTau = 
        mesh().template lookupObject<volScalarField>("rDeltaT");

    fvm.diag() = (rDeltaTau.primitiveField() + coefft*rDeltaT)*mesh().V();

    if (mesh().moving())
    {
        fvm.source() =
        (
            rDeltaTau*vf.primitiveField()*mesh().V()
          + coefft0*rDeltaT*vf.oldTime().primitiveField()*mesh().V0()
          - coefft00*rDeltaT*vf.oldTime().oldTime().primitiveField()
           *mesh().V00()
        );
    }
    else
    {
        fvm.source() = mesh().V()*
        (
            rDeltaTau*vf.primitiveField()
          + coefft0*rDeltaT*vf.oldTime().primitiveField()
          - coefft00*rDeltaT*vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
backwardDualTimeDdtScheme<Type>::fvmDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    const volScalarField& rDeltaTau = 
        mesh().template lookupObject<volScalarField>("rDeltaT");

    fvm.diag() = 
        (rDeltaTau.primitiveField() + coefft*rDeltaT)
       *rho.value()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rho.value()*
        (
            rDeltaTau*vf.primitiveField()*mesh().V()
          + coefft0*rDeltaT*vf.oldTime().primitiveField()*mesh().V0()
          - coefft00*rDeltaT*vf.oldTime().oldTime().primitiveField()
           *mesh().V00()
        );
    }
    else
    {
        fvm.source() = mesh().V()*rho.value()*
        (
            rDeltaTau*vf.primitiveField()
          + coefft0*rDeltaT*vf.oldTime().primitiveField()
          - coefft00*rDeltaT*vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
backwardDualTimeDdtScheme<Type>::fvmDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    const volScalarField& rDeltaTau = 
        mesh().template lookupObject<volScalarField>("rDeltaT");

    fvm.diag() =
        (rDeltaTau.primitiveField() + coefft*rDeltaT)
       *rho.primitiveField()*mesh().V();
 
    if (mesh().moving())
    {
        fvm.source() = 
        (
            rDeltaTau
           *rho.prevIter().primitiveField()
           *vf.primitiveField()*mesh().V()

          + coefft0*rDeltaT
           *rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()*mesh().V0()

          - coefft00*rDeltaT
           *rho.oldTime().oldTime().primitiveField()
           *vf.oldTime().oldTime().primitiveField()*mesh().V00()
        );
    }
    else
    {
        fvm.source() = mesh().V()*
        (
            rDeltaTau
           *rho.prevIter().primitiveField()
           *vf.primitiveField()

          + coefft0*rDeltaT
           *rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()

          - coefft00*rDeltaT
           *rho.oldTime().oldTime().primitiveField()
           *vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
backwardDualTimeDdtScheme<Type>::fvmDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    tmp<fvMatrix<Type>> tfvm
    (
        new fvMatrix<Type>
        (
            vf,
            alpha.dimensions()*rho.dimensions()*vf.dimensions()*dimVol/dimTime
        )
    );
    fvMatrix<Type>& fvm = tfvm.ref();

    scalar rDeltaT = 1.0/deltaT_();

    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    scalar coefft   = 1 + deltaT/(deltaT + deltaT0);
    scalar coefft00 = deltaT*deltaT/(deltaT0*(deltaT + deltaT0));
    scalar coefft0  = coefft + coefft00;

    const volScalarField& rDeltaTau = 
        mesh().template lookupObject<volScalarField>("rDeltaT");

    fvm.diag() =
        (rDeltaTau.primitiveField() + coefft*rDeltaT)
       *alpha.primitiveField()*rho.primitiveField()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = 
        (
            rDeltaTau
           *alpha.prevIter().primitiveField()
           *rho.prevIter().primitiveField()
           *vf.primitiveField()*mesh().V()

          + coefft0*rDeltaT
           *alpha.oldTime().primitiveField()
           *rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()*mesh().V0()

          - coefft00*rDeltaT
           *alpha.oldTime().oldTime().primitiveField()
           *rho.oldTime().oldTime().primitiveField()
           *vf.oldTime().oldTime().primitiveField()*mesh().V00()
        );
    }
    else
    {
        fvm.source() = mesh().V()*
        (
            rDeltaTau
           *alpha.prevIter().primitiveField()
           *rho.prevIter().primitiveField()
           *vf.primitiveField()

          + coefft0*rDeltaT
           *alpha.oldTime().primitiveField()
           *rho.oldTime().primitiveField()
           *vf.oldTime().primitiveField()

          - coefft00*rDeltaT
           *alpha.oldTime().oldTime().primitiveField()
           *rho.oldTime().oldTime().primitiveField()
           *vf.oldTime().oldTime().primitiveField()
        );
    }

    return tfvm;
}


template<class Type>
tmp<typename backwardDualTimeDdtScheme<Type>::fluxFieldType>
backwardDualTimeDdtScheme<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    NotImplemented;
    return fluxFieldType::null();
}


template<class Type>
tmp<typename backwardDualTimeDdtScheme<Type>::fluxFieldType>
backwardDualTimeDdtScheme<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    NotImplemented;
    return fluxFieldType::null();
}


template<class Type>
tmp<typename backwardDualTimeDdtScheme<Type>::fluxFieldType>
backwardDualTimeDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    NotImplemented;
    return fluxFieldType::null();
}


template<class Type>
tmp<typename backwardDualTimeDdtScheme<Type>::fluxFieldType>
backwardDualTimeDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    
    NotImplemented;
    return fluxFieldType::null();
}


template<class Type>
tmp<surfaceScalarField> backwardDualTimeDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    scalar deltaT = deltaT_();
    scalar deltaT0 = deltaT0_(vf);

    // Coefficient for t-3/2 (between times 0 and 00)
    scalar coefft0_00 = deltaT/(deltaT + deltaT0);

    // Coefficient for t-1/2 (between times n and 0)
    scalar coefftn_0 = 1 + coefft0_00;

    return tmp<surfaceScalarField>
    (
        new surfaceScalarField
        (
            IOobject
            (
                mesh().phi().name(),
                mesh().time().timeName(),
                mesh().thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            coefftn_0*mesh().phi() - coefft0_00*mesh().phi().oldTime()
        )
    );
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
