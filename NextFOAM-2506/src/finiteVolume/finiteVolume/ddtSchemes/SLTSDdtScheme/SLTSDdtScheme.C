/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2018 OpenFOAM Foundation
    Copyright (C) 2023 OpenCFD Ltd.
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

#include "SLTSDdtScheme.H"
#include "surfaceInterpolate.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
void SLTSDdtScheme<Type>::relaxedDiag
(
    scalarField& rD,
    const surfaceScalarField& phi
) const
{
    const labelUList& owner = mesh().owner();
    const labelUList& neighbour = mesh().neighbour();
    scalarField diag(rD.size(), Zero);

    forAll(owner, facei)
    {
        if (phi[facei] > 0.0)
        {
            diag[owner[facei]] += phi[facei];
            rD[neighbour[facei]] += phi[facei];
        }
        else
        {
            diag[neighbour[facei]] -= phi[facei];
            rD[owner[facei]] -= phi[facei];
        }
    }

    forAll(phi.boundaryField(), patchi)
    {
        const fvsPatchScalarField& pphi = phi.boundaryField()[patchi];
        const labelUList& faceCells = pphi.patch().patch().faceCells();

        forAll(pphi, patchFacei)
        {
            if (pphi[patchFacei] > 0.0)
            {
                diag[faceCells[patchFacei]] += pphi[patchFacei];
            }
            else
            {
                rD[faceCells[patchFacei]] -= pphi[patchFacei];
            }
        }
    }

    rD += (1.0/alpha_ - 2.0)*diag;
}


template<class Type>
tmp<volScalarField> SLTSDdtScheme<Type>::SLrDeltaT() const
{
    const surfaceScalarField& phi =
        mesh().objectRegistry::template
            lookupObject<surfaceScalarField>(phiName_);

    const dimensionedScalar& deltaT = mesh().time().deltaT();

    tmp<volScalarField> trDeltaT
    (
        new volScalarField
        (
            IOobject
            (
                "rDeltaT",
                phi.instance(),
                mesh().thisDb()
            ),
            mesh(),
            dimensionedScalar(dimless/dimTime, Zero),
            fvPatchFieldBase::extrapolatedCalculatedType()
        )
    );

    volScalarField& rDeltaT = trDeltaT.ref();

    relaxedDiag(rDeltaT, phi);

    if (phi.dimensions() == dimensionSet(0, 3, -1, 0, 0))
    {
        rDeltaT.primitiveFieldRef() = max
        (
            rDeltaT.primitiveField()/mesh().V(),
            scalar(1)/deltaT.value()
        );
    }
    else if (phi.dimensions() == dimensionSet(1, 0, -1, 0, 0))
    {
        const volScalarField& rho =
            mesh().objectRegistry::template lookupObject<volScalarField>
            (
                rhoName_
            ).oldTime();

        rDeltaT.primitiveFieldRef() = max
        (
            rDeltaT.primitiveField()/(rho.primitiveField()*mesh().V()),
            scalar(1)/deltaT.value()
        );
    }
    else
    {
        FatalErrorInFunction
            << "Incorrect dimensions of phi: " << phi.dimensions()
            << abort(FatalError);
    }

    rDeltaT.correctBoundaryConditions();

    return trDeltaT;
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
SLTSDdtScheme<Type>::fvcDdt
(
    const dimensioned<Type>& dt
)
{
    const volScalarField rDeltaT(SLrDeltaT());

    IOobject ddtIOobject
    (
        "ddt("+dt.name()+')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

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

        tdtdt.ref().primitiveFieldRef() =
            rDeltaT.primitiveField()*dt.value()*(1.0 - mesh().V0()/mesh().V());

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
                mesh(),
                dimensioned<Type>(dt.dimensions()/dimTime, Zero),
                fvPatchFieldBase::calculatedType()
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
SLTSDdtScheme<Type>::fvcDdt
(
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField rDeltaT(SLrDeltaT());

    IOobject ddtIOobject
    (
        "ddt("+vf.name()+')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    if (mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*vf.dimensions(),
                rDeltaT.primitiveField()*
                (
                    vf.primitiveField()
                  - vf.oldTime().primitiveField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.boundaryField()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
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
                rDeltaT*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
SLTSDdtScheme<Type>::fvcDdt
(
    const dimensionedScalar& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField rDeltaT(SLrDeltaT());

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    if (mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.primitiveField()*rho.value()*
                (
                    vf.primitiveField()
                  - vf.oldTime().primitiveField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.boundaryField()*rho.value()*
                (
                    vf.boundaryField() - vf.oldTime().boundaryField()
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
                rDeltaT*rho*(vf - vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
SLTSDdtScheme<Type>::fvcDdt
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField rDeltaT(SLrDeltaT());

    IOobject ddtIOobject
    (
        "ddt("+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

    if (mesh().moving())
    {
        tmp<GeometricField<Type, fvPatchField, volMesh>> tdtdt
        (
            new GeometricField<Type, fvPatchField, volMesh>
            (
                ddtIOobject,
                mesh(),
                rDeltaT.dimensions()*rho.dimensions()*vf.dimensions(),
                rDeltaT.primitiveField()*
                (
                    rho.primitiveField()*vf.primitiveField()
                  - rho.oldTime().primitiveField()
                   *vf.oldTime().primitiveField()*mesh().V0()/mesh().V()
                ),
                rDeltaT.boundaryField()*
                (
                    rho.boundaryField()*vf.boundaryField()
                  - rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
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
                rDeltaT*(rho*vf - rho.oldTime()*vf.oldTime())
            )
        );
    }
}


template<class Type>
tmp<GeometricField<Type, fvPatchField, volMesh>>
SLTSDdtScheme<Type>::fvcDdt
(
    const volScalarField& alpha,
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& vf
)
{
    const volScalarField rDeltaT(SLrDeltaT());

    IOobject ddtIOobject
    (
        "ddt("+alpha.name()+','+rho.name()+','+vf.name()+')',
        mesh().time().timeName(),
        mesh().thisDb()
    );

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
                rDeltaT.primitiveField()*
                (
                    alpha.primitiveField()
                   *rho.primitiveField()
                   *vf.primitiveField()

                  - alpha.oldTime().primitiveField()
                   *rho.oldTime().primitiveField()
                   *vf.oldTime().primitiveField()*mesh().Vsc0()/mesh().Vsc()
                ),
                rDeltaT.boundaryField()*
                (
                    alpha.boundaryField()
                   *rho.boundaryField()
                   *vf.boundaryField()

                  - alpha.oldTime().boundaryField()
                   *rho.oldTime().boundaryField()
                   *vf.oldTime().boundaryField()
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
                rDeltaT
               *(
                   alpha*rho*vf
                 - alpha.oldTime()*rho.oldTime()*vf.oldTime()
                )
            )
        );
    }
}


template<class Type>
tmp<fvMatrix<Type>>
SLTSDdtScheme<Type>::fvmDdt
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

    scalarField rDeltaT(SLrDeltaT()().primitiveField());

    fvm.diag() = rDeltaT*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT*vf.oldTime().primitiveField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT*vf.oldTime().primitiveField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
SLTSDdtScheme<Type>::fvmDdt
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

    scalarField rDeltaT(SLrDeltaT()().primitiveField());

    fvm.diag() = rDeltaT*rho.value()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().primitiveField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.value()*vf.oldTime().primitiveField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
SLTSDdtScheme<Type>::fvmDdt
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

    scalarField rDeltaT(SLrDeltaT()().primitiveField());

    fvm.diag() = rDeltaT*rho.primitiveField()*mesh().V();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().V0();
    }
    else
    {
        fvm.source() = rDeltaT
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().V();
    }

    return tfvm;
}


template<class Type>
tmp<fvMatrix<Type>>
SLTSDdtScheme<Type>::fvmDdt
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

    scalarField rDeltaT(SLrDeltaT()().primitiveField());

    fvm.diag() =
        rDeltaT*alpha.primitiveField()*rho.primitiveField()*mesh().Vsc();

    if (mesh().moving())
    {
        fvm.source() = rDeltaT
            *alpha.oldTime().primitiveField()
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().Vsc0();
    }
    else
    {
        fvm.source() = rDeltaT
            *alpha.oldTime().primitiveField()
            *rho.oldTime().primitiveField()
            *vf.oldTime().primitiveField()*mesh().Vsc();
    }

    return tfvm;
}


template<class Type>
tmp<typename SLTSDdtScheme<Type>::fluxFieldType>
SLTSDdtScheme<Type>::fvcDdtUfCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    const surfaceScalarField rDeltaT(fvc::interpolate(SLrDeltaT()));

    fluxFieldType phiUf0(mesh().Sf() & Uf.oldTime());
    fluxFieldType phiCorr
    (
        phiUf0 - fvc::dotInterpolate(mesh().Sf(), U.oldTime())
    );

    return tmp<fluxFieldType>
    (
        new fluxFieldType
        (
            IOobject
            (
                "ddtCorr(" + U.name() + ',' + Uf.name() + ')',
                mesh().time().timeName(),
                mesh().thisDb()
            ),
            this->fvcDdtPhiCoeff(U.oldTime(), phiUf0, phiCorr)
           *rDeltaT*phiCorr
        )
    );
}


template<class Type>
tmp<typename SLTSDdtScheme<Type>::fluxFieldType>
SLTSDdtScheme<Type>::fvcDdtPhiCorr
(
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    if (mesh().objectRegistry::template foundObject<volScalarField>("rUAp"))
    {
        volScalarField rDeltaT(SLrDeltaT());

        tmp<volScalarField> trAU
        (
            new volScalarField
            (
                IOobject
                (
                    "rUAp",
                    mesh().time().timeName(),
                    mesh().thisDb()
                ),
                mesh().objectRegistry::template lookupObject<volScalarField>
                (
                    "rUAp"
                )
            )
        );

        const volScalarField& rAU = trAU();

        volScalarField dAU("dAU", rAU*rDeltaT);
        surfaceScalarField dAUf("dAUf", momentumInterpolate(dAU));

        tmp<fluxFieldType> toldUPhi
        (
            momentumInterpolate(dAU*U.oldTime()) & mesh().Sf()
        );

        tmp<fluxFieldType> toldPhi
        (
            dAUf*phi.oldTime()
        );

        tmp<fluxFieldType> tPhiCorr
        (
            new fluxFieldType
            (
                IOobject
                (
                    "ddtCorr(" + U.name() + ',' + phi.name() + ')',
                    mesh().time().timeName(),
                    mesh().thisDb()
                ),
                mesh(),
                dimensioned<typename flux<Type>::type>
                (
                    "zero",
                    phi.dimensions(),
                    pTraits<typename flux<Type>::type>::zero
                )
            )
        );

        Field<typename flux<Type>::type>& phiCorr
        (
            tPhiCorr.ref().primitiveFieldRef()
        );

        phiCorr = toldPhi().primitiveField() - toldUPhi().primitiveField();

        typename GeometricField
        <
            typename flux<Type>::type, fvsPatchField, surfaceMesh
        >::Boundary& bPhiCorr = tPhiCorr.ref().boundaryFieldRef();

        const volScalarField& p
        (
            mesh().objectRegistry::template foundObject<volScalarField>
            (
                "p_rgh"
            )
          ? mesh().objectRegistry::template lookupObject<volScalarField>
            (
                "p_rgh"
            )
          : mesh().objectRegistry::template lookupObject<volScalarField>("p")
        );

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
                    toldPhi().boundaryField()[patchi]
                  - toldUPhi().boundaryField()[patchi];
            }
        }

        return tPhiCorr;
    } // by Gill
    else
    {
        const surfaceScalarField rDeltaT(fvc::interpolate(SLrDeltaT()));

        fluxFieldType phiCorr
        (
            phi.oldTime() - fvc::dotInterpolate(mesh().Sf(), U.oldTime())
        );

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                IOobject
                (
                    "ddtCorr(" + U.name() + ',' + phi.name() + ')',
                    mesh().time().timeName(),
                    mesh().thisDb()
                ),
                this->fvcDdtPhiCoeff(U.oldTime(), phi.oldTime(), phiCorr)
               *rDeltaT*phiCorr
            )
        );
    }
}


template<class Type>
tmp<typename SLTSDdtScheme<Type>::fluxFieldType>
SLTSDdtScheme<Type>::fvcDdtUfCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const GeometricField<Type, fvsPatchField, surfaceMesh>& Uf
)
{
    const surfaceScalarField rDeltaT(fvc::interpolate(SLrDeltaT()));

    if
    (
        U.dimensions() == dimVelocity
     && Uf.dimensions() == dimDensity*dimVelocity
    )
    {
        GeometricField<Type, fvPatchField, volMesh> rhoU0
        (
            rho.oldTime()*U.oldTime()
        );

        fluxFieldType phiUf0(mesh().Sf() & Uf.oldTime());
        fluxFieldType phiCorr(phiUf0 - fvc::dotInterpolate(mesh().Sf(), rhoU0));

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                IOobject
                (
                    "ddtCorr("
                  + rho.name() + ',' + U.name() + ',' + Uf.name() + ')',
                    mesh().time().timeName(),
                    mesh().thisDb()
                ),
                this->fvcDdtPhiCoeff(rhoU0, phiUf0, phiCorr, rho.oldTime())
               *rDeltaT*phiCorr
            )
        );
    }
    else if
    (
        U.dimensions() == dimDensity*dimVelocity
     && Uf.dimensions() == dimDensity*dimVelocity
    )
    {
        fluxFieldType phiUf0(mesh().Sf() & Uf.oldTime());
        fluxFieldType phiCorr
        (
            phiUf0 - fvc::dotInterpolate(mesh().Sf(), U.oldTime())
        );

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                IOobject
                (
                    "ddtCorr("
                  + rho.name() + ',' + U.name() + ',' + Uf.name() + ')',
                    mesh().time().timeName(),
                    mesh().thisDb()
                ),
                this->fvcDdtPhiCoeff
                (
                    U.oldTime(),
                    phiUf0,
                    phiCorr,
                    rho.oldTime()
                )*rDeltaT*phiCorr
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of Uf are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<typename SLTSDdtScheme<Type>::fluxFieldType>
SLTSDdtScheme<Type>::fvcDdtPhiCorr
(
    const volScalarField& rho,
    const GeometricField<Type, fvPatchField, volMesh>& U,
    const fluxFieldType& phi
)
{
    if
    (
        U.dimensions() == dimVelocity
     && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        if 
        (
            mesh().objectRegistry::template foundObject<volScalarField>("rUAp")
        )
        {
            volScalarField rDeltaT(SLrDeltaT());

            tmp<volScalarField> trAU
            (
                new volScalarField
                (
                    IOobject
                    (
                        "rUAp",
                        mesh().time().timeName(),
                        mesh().thisDb()
                    ),
                    mesh().objectRegistry::template lookupObject
                    <
                        volScalarField
                    >("rUAp")
                )
            );

            const volScalarField& rAU = trAU();

            volScalarField dAU("dAU", rAU*rDeltaT);
            surfaceScalarField dAUf("dAUf", momentumInterpolate(dAU));

            const surfaceScalarField& rhof
            (
                mesh().objectRegistry::template lookupObject
                <
                    surfaceScalarField
                >("faceRho")
            );

            tmp<fluxFieldType> toldUPhi
            (
                rhof*(momentumInterpolate(dAU*U.oldTime()) & mesh().Sf())
            );

            tmp<fluxFieldType> toldPhi(dAUf*phi.oldTime());

            tmp<fluxFieldType> tPhiCorr
            (
                new fluxFieldType
                (
                    IOobject
                    (
                        "ddtCorr("
                      + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
                        mesh().time().timeName(),
                        mesh().thisDb()
                    ),
                    mesh(),
                    dimensioned<typename flux<Type>::type>
                    (
                        "zero",
                        phi.dimensions(),
                        pTraits<typename flux<Type>::type>::zero
                    )
                )
            );

            Field<typename flux<Type>::type>& phiCorr
            (
                tPhiCorr.ref().primitiveFieldRef()
            );

            phiCorr = toldPhi().primitiveField() - toldUPhi().primitiveField();

            typename GeometricField
            <
                typename flux<Type>::type, fvsPatchField, surfaceMesh
            >::Boundary& bPhiCorr = tPhiCorr.ref().boundaryFieldRef();

            const volScalarField& p
            (
                mesh().objectRegistry::template foundObject<volScalarField>
                (
                    "p_rgh"
                )
              ? mesh().objectRegistry::template lookupObject<volScalarField>
                (
                    "p_rgh"
                )
              : mesh().objectRegistry::template lookupObject<volScalarField>
                (
                    "p"
                )
            );

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
                    (
                       toldPhi().boundaryField()[patchi]
                     - toldUPhi().boundaryField()[patchi]
                    );
                }
            }

            return tPhiCorr;
        } // by Gill
        else
        {
            dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

            GeometricField<Type, fvPatchField, volMesh> rhoU0
            (
                rho.oldTime()*U.oldTime()
            );

            fluxFieldType phiCorr
            (
                phi.oldTime() - fvc::dotInterpolate(mesh().Sf(), rhoU0)
            );

            return tmp<fluxFieldType>
            (
                new fluxFieldType
                (
                    IOobject
                    (
                        "ddtCorr("
                      + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
                        mesh().time().timeName(),
                        mesh().thisDb()
                    ),
                    this->fvcDdtPhiCoeff
                    (
                        rhoU0,
                        phi.oldTime(),
                        phiCorr,
                        rho.oldTime()
                    )*rDeltaT*phiCorr
                )
            );
        }
    }
    else if
    (
        U.dimensions() == rho.dimensions()*dimVelocity
     && phi.dimensions() == rho.dimensions()*dimVelocity*dimArea
    )
    {
        dimensionedScalar rDeltaT = 1.0/mesh().time().deltaT();

        fluxFieldType phiCorr
        (
            phi.oldTime() - fvc::dotInterpolate(mesh().Sf(), U.oldTime())
        );

        return tmp<fluxFieldType>
        (
            new fluxFieldType
            (
                IOobject
                (
                    "ddtCorr("
                  + rho.name() + ',' + U.name() + ',' + phi.name() + ')',
                    mesh().time().timeName(),
                    mesh().thisDb()
                ),
                this->fvcDdtPhiCoeff
                (
                    U.oldTime(),
                    phi.oldTime(),
                    phiCorr,
                    rho.oldTime()
                )*rDeltaT*phiCorr
            )
        );
    }
    else
    {
        FatalErrorInFunction
            << "dimensions of phi are not correct"
            << abort(FatalError);

        return fluxFieldType::null();
    }
}


template<class Type>
tmp<surfaceScalarField> SLTSDdtScheme<Type>::meshPhi
(
    const GeometricField<Type, fvPatchField, volMesh>&
)
{
    tmp<surfaceScalarField> tmeshPhi
    (
        new surfaceScalarField
        (
            IOobject
            (
                "meshPhi",
                mesh().time().timeName(),
                mesh().thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                IOobject::NO_REGISTER
            ),
            mesh(),
            dimensionedScalar(dimVolume/dimTime, Zero)
        )
    );

    tmeshPhi.ref().setOriented();

    return tmeshPhi;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
