/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "realizableKEtwoLayer.H"
#include "fvOptions.H"
#include "bound.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
tmp<volScalarField> realizableKEtwoLayer<BasicTurbulenceModel>::rCmu
(
    const volTensorField& gradU,
    const volScalarField& S2,
    const volScalarField& magS
)
{
    tmp<volSymmTensorField> tS = dev(symm(gradU));
    const volSymmTensorField& S = tS();

    tmp<volScalarField> magSS2(magS*S2);
    magSS2.ref().clamp_min(VSMALL);
    volScalarField W((2*sqrt(2.0))*((S&S)&&S)/magSS2);

    tS.clear();

    volScalarField phis
    (
        (1.0/3.0)*acos(clamp(sqrt(6.0)*W, scalarMinMax(-1, 1)))
    );
    volScalarField As(sqrt(6.0)*cos(phis));
    volScalarField Us(sqrt(S2/2.0 + magSqr(skew(gradU))));

    epsilon_.clamp_min(VSMALL);
    tmp<volScalarField> rrCmu(A0_ + As*Us*k_/epsilon_);
    rrCmu.ref().clamp_min(VSMALL);

    return 1.0/rrCmu;
}


template<class BasicTurbulenceModel>
void realizableKEtwoLayer<BasicTurbulenceModel>::correctNut
(
    const volTensorField& gradU,
    const volScalarField& S2,
    const volScalarField& magS
)
{
    this->nut_ =
        max
        (
            min
            (
                rCmu(gradU, S2, magS)*sqr(k_)/epsilon_,
                this->viscosityRatioMax_*this->nu()
            ),
            this->viscosityRatioMin_*this->nu()
        ); 

    dimensionedScalar A = deltaRey_/tanh(0.98);

    Rey_ = y_()*sqrt(this->k_())/this->nu()();
    lambda_ = 0.5*(1.0 + tanh((Rey_ - ReyStar_)/A));
    //lambda_.clamp_min(VSMALL);

    volScalarField::Internal lmu(y_()*kappa_/Cmu75_*(1 - exp(-Rey_/Anu_)));
    volScalarField::Internal nut_tl(Cmu_*lmu*sqrt(k_()));
    const_cast<volScalarField::Internal&>(this->nut_()) = 
        lambda_*this->nut_() + (1 - lambda_)*nut_tl;

    this->nut_.correctBoundaryConditions();
    fv::options::New(this->mesh_).correct(this->nut_);

    BasicTurbulenceModel::correctNut();
}


template<class BasicTurbulenceModel>
void realizableKEtwoLayer<BasicTurbulenceModel>::correctNut()
{
    tmp<volTensorField> tgradU = fvc::grad(this->U_);
    volScalarField S2(2*magSqr(dev(symm(tgradU()))));
    volScalarField magS(sqrt(S2));
    correctNut(tgradU(), S2, magS);
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> realizableKEtwoLayer<BasicTurbulenceModel>::kSource() const
{
    // for buoyancy
    if (this->mesh_.template foundObject<uniformDimensionedVectorField>("g"))
    {
        const uniformDimensionedVectorField& g =
            this->mesh_.template
            lookupObject<uniformDimensionedVectorField>("g");

        if (mag(g.value()) > SMALL && this->rho_.dimensions() == dimDensity)
        {
            return
                -fvm::SuSp
                (
                    this->nut_
                    *(
                        g
                      & fvc::grad
                        (
                            reinterpret_cast<const volScalarField&>(this->rho_)
                        )
                     )/Prt_/this->k_,
                    this->k_
                );
        }
        else
        {
            return tmp<fvScalarMatrix>
            (
                new fvScalarMatrix
                (
                    k_,
                    dimVolume*this->rho_.dimensions()*k_.dimensions()
                    /dimTime
                )
            );
        }
    }
    else
    {
        return tmp<fvScalarMatrix>    
        (
            new fvScalarMatrix
            (
                k_,              
                dimVolume*this->rho_.dimensions()*k_.dimensions()
                /dimTime
            )    
        );
    }
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix> realizableKEtwoLayer<BasicTurbulenceModel>::epsilonSource() const
{
    return tmp<fvScalarMatrix>
    (
        new fvScalarMatrix
        (
            epsilon_,
            dimVolume*this->rho_.dimensions()*epsilon_.dimensions()
            /dimTime
        )
    );
}


template<class BasicTurbulenceModel>
void realizableKEtwoLayer<BasicTurbulenceModel>::blendAndRelax
(
    fvScalarMatrix& eqn,
    const volScalarField::Internal& lambda,
    const volScalarField::Internal& es
)
{
    scalar alpha(1.0);

    word name = epsilon_.select
    (
        epsilon_.mesh().data().isFinalIteration()
    ); 

    if (epsilon_.mesh().relaxEquation(name))
    {
        alpha = epsilon_.mesh().equationRelaxationFactor(name);
    }

    if (alpha <= 0)
    {
        return;
    }

    DebugInFunction
        << "Relaxing " << epsilon_.name() << " by " << alpha << endl;

    // Get unrelaxed central coefficient
    tmp<scalarField> scaledAp = eqn.D()*(1 - lambda)/lambda;

    // Calculate the sum-mag off-diagonal from the interior faces
    scalarField sumOff(scaledAp().size(), Zero);
    eqn.sumMagOffDiag(sumOff);

    scalarField& S = eqn.source();
    scalarField& D = eqn.diag();

    // Store the current unrelaxed diagonal for use in updating the source
    scalarField D0(D);

    // Handle the boundary contributions to the diagonal
    forAll(epsilon_.boundaryField(), patchi)
    {
        const fvPatchScalarField& ptf = epsilon_.boundaryField()[patchi];

        if (ptf.size())
        {
            const labelUList& pa = eqn.lduAddr().patchAddr(patchi);
            const scalarField& iCoeffs = eqn.internalCoeffs()[patchi];

            if (ptf.coupled())
            {
                const scalarField& pCoeffs = eqn.boundaryCoeffs()[patchi];

  
                // For coupled boundaries add the diagonal and
                // off-diagonal contributions              
                forAll(pa, face)
                {
                    D[pa[face]] += iCoeffs[face];
                    sumOff[pa[face]] += mag(pCoeffs[face]);
                }               
            }
            else
            {
                // For non-coupled boundaries add the magnitude of diagonal
                // contribution to ensure stability               
                forAll(pa, face)
                {
                    D[pa[face]] += mag(iCoeffs[face]);
                }               
            }
        }
    }

    // Ensure the matrix is diagonally dominant...
    // Assumes that the central coefficient is positive and ensures it is
    forAll(D, celli)
    {
        D[celli] = max(mag(D[celli]), sumOff[celli]);
    }

    // ... then relax
    D /= (alpha*lambda);

    // Now remove the diagonal contribution from the boundaries
    forAll(epsilon_.boundaryField(), patchi)
    {
        const fvPatchScalarField& ptf = epsilon_.boundaryField()[patchi];

        if (ptf.size())
        {
            const labelUList& pa = eqn.lduAddr().patchAddr(patchi);
            const scalarField& iCoeffs = eqn.internalCoeffs()[patchi];

            forAll(pa, face)
            {
                D[pa[face]] -= iCoeffs[face];
            }
        }
    }

    // Finally add the relaxation and blending contribution to the source.
    S += (D - D0 - scaledAp())*epsilon_.primitiveField();
    S += scaledAp*es;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
realizableKEtwoLayer<BasicTurbulenceModel>::realizableKEtwoLayer
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<RASModel<BasicTurbulenceModel>>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),
    A0_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "A0",
            this->coeffDict_,
            4.04
        )
    ),
    C2_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "C2",
            this->coeffDict_,
            1.9
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmak",
            this->coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "sigmaEps",
            this->coeffDict_,
            1.2
        )
    ),
    Cmu_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Cmu",
            this->coeffDict_,
            0.09
        )
    ),
    Cmu75_(pow(Cmu_, 0.75)),
    kappa_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "kappa",
            this->coeffDict_,
            0.41
        )
    ),
    ReyStar_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "ReyStar",
            this->coeffDict_,
            60
        )
    ),
    deltaRey_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "deltaRey",
            this->coeffDict_,
            10
        )
    ),
    Anu_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Anu",
            this->coeffDict_,
            70
        )
    ),
    Prt_
    (
        dimensioned<scalar>::getOrAddToDict
        (
            "Prt",
            this->coeffDict_,
            0.85
        )
    ),

    k_
    (
        IOobject
        (
            IOobject::groupName("k", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    epsilon_
    (
        IOobject
        (
            IOobject::groupName("epsilon", U.group()),
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        this->mesh_
    ),
    rCmu_
    (
        IOobject
        (
            "rCmu",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        Cmu_
    ),
    lambda_
    (
        IOobject
        (
            "lambda",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("one", dimless, 1.0)
    ),
    Rey_
    (
        IOobject
        (
            "Rey",
            this->runTime_.timeName(),
            this->mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        this->mesh_,
        dimensionedScalar("one", dimless, 1.0)
    ),
    y_(wallDist::New(this->mesh_).y())
{
    bound(k_, this->kMin_);
    bound(epsilon_, this->epsilonMin_);

    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool realizableKEtwoLayer<BasicTurbulenceModel>::read()
{
    if (eddyViscosity<RASModel<BasicTurbulenceModel>>::read())
    {
        A0_.readIfPresent(this->coeffDict());
        C2_.readIfPresent(this->coeffDict());
        sigmak_.readIfPresent(this->coeffDict());
        sigmaEps_.readIfPresent(this->coeffDict());
        Cmu_.readIfPresent(this->coeffDict());
        Cmu75_ = pow(Cmu_, 0.75),
        kappa_.readIfPresent(this->coeffDict());
        ReyStar_.readIfPresent(this->coeffDict());
        deltaRey_.readIfPresent(this->coeffDict());
        Anu_.readIfPresent(this->coeffDict());

        return true;
    }

    return false;
}


template<class BasicTurbulenceModel>
void realizableKEtwoLayer<BasicTurbulenceModel>::correct()
{
    if (!this->turbulence_)
    {
        return;
    }

    // Local references
    const alphaField& alpha = this->alpha_;
    const rhoField& rho = this->rho_;
    const surfaceScalarField& alphaRhoPhi = this->alphaRhoPhi_;
    const volVectorField& U = this->U_;
    const volScalarField& nut = this->nut_;

    fv::options& fvOptions(fv::options::New(this->mesh_));

    eddyViscosity<RASModel<BasicTurbulenceModel>>::correct();

    tmp<volScalarField::Internal> divU = turbulenceModel::divU();

    tmp<volTensorField> tgradU = fvc::grad(U);
    const volScalarField::Internal GbyNu
    (
        this->type() + ":GbyNu",
        tgradU().v() && dev(twoSymm(tgradU().v()))
    );
    //const volScalarField::Internal G(this->GName(), nut()*GbyNu);
    const volScalarField::Internal G(this->GName(), this->nuEff()()*GbyNu);

    volScalarField S2("S2", 2*magSqr(dev(symm(tgradU()))));
    volScalarField magS(sqrt(S2));

    volScalarField eta(magS*k_/epsilon_);

    tmp<volScalarField> etaPlus5(scalar(5) + eta);
    etaPlus5.ref().clamp_min(VSMALL);

    volScalarField C1(max(eta/etaPlus5, scalar(0.43)));

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();
    // Push any changed cell values to coupled neighbours
    epsilon_.boundaryFieldRef().template evaluateCoupled<coupledFvPatch>();

    rCmu_ = rCmu(tgradU(), S2, magS);

    // SAF: limiting thermo->nu(). If psiThermo is used rho might be < 0
    // temporarily when p < 0 then nu < 0 which needs limiting
    const volScalarField::Internal rootNuLimited
    (
        sqrt
        (
            max
            (
                this->nu()(),
                dimensionedScalar(this->nu()().dimensions(), Zero)
            )
        )
    );

    const volScalarField::Internal rootNutByCmu(sqrt(nut()/rCmu_()));

    // Dissipation equation - modified source term linearization
    const volScalarField::Internal epsilonSourceCoeff
    (
        0.5*alpha()*rho()*sqrt(epsilon_())
        *(
            C2_/(rootNutByCmu + rootNuLimited) - C1*eta()/rootNutByCmu
         )
    );

    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(alpha, rho, epsilon_)
      + fvm::div(alphaRhoPhi, epsilon_)
      - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
     ==
        epsilonSourceCoeff*epsilon_()
      - fvm::SuSp(3.0*epsilonSourceCoeff, epsilon_)
      + epsilonSource()
      + fvOptions(alpha, rho, epsilon_)
    );

    // Blend with one-equation turbulence model near the wall
    {
        dimensionedScalar Aeps = 2.0*kappa_/Cmu75_;
        volScalarField::Internal leps
        (
            y_()*kappa_/Cmu75_*(1 - exp(-Rey_/Aeps))
        );
        volScalarField::Internal es(pow(k_(), 1.5)/leps);

        blendAndRelax(epsEqn.ref(), lambda_, es);
    }

    fvOptions.constrain(epsEqn.ref());
    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());
    solve(epsEqn);
    fvOptions.correct(epsilon_);
    bound(epsilon_, this->epsilonMin_);

    // Turbulent kinetic energy equation - modified source term linearization
    const volScalarField::Internal kSourceCoeff
    (
        alpha()*rho()*rCmu_()*(k_()/nut())
    );

    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(alpha, rho, k_)
      + fvm::div(alphaRhoPhi, k_)
      - fvm::laplacian(alpha*rho*DkEff(), k_)
     ==
        alpha*rho*G
      - fvm::SuSp((2.0/3.0)*alpha()*rho()*divU, k_)
      + kSourceCoeff*k_()
      - fvm::SuSp(2.0*kSourceCoeff, k_)
      + kSource()
      + fvOptions(alpha, rho, k_)
    );

    kEqn.ref().relax();
    fvOptions.constrain(kEqn.ref());
    solve(kEqn);
    fvOptions.correct(k_);
    bound(k_, this->kMin_);

    correctNut(tgradU(), S2, magS);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
