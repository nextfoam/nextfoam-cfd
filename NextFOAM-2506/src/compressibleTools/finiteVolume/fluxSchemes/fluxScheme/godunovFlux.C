/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
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

#include "godunovFlux.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace compressible
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(godunovFlux, 0);
defineRunTimeSelectionTable(godunovFlux, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

godunovFlux::godunovFlux
(
	const word& type, 
    const volScalarField& rho,
    const volVectorField& U,
    const fluidThermo& thermoPhysicalModel,
    const turbulenceModel& turbulenceModel
)
:
	dictionary(U.mesh().solutionDict().subDict("Riemann")),
    runTime_(U.time()),
    mesh_(U.mesh()),

    rho_(rho),
    U_(U),

    p_(thermoPhysicalModel.p()),
    he_(thermoPhysicalModel.he()),

    thermoPhysicalModel_(thermoPhysicalModel),
    turbulenceModel_(turbulenceModel),

    // fluxes in mass conservation equation: \varrho \vec{u}
    // only initialization!
    rhoFlux_
    (
        IOobject
        (
            "rhoFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        (linearInterpolate(rho_*U_) & mesh_.Sf())
    ),
    // fluxes in momentum equation: \varrho \vec{u} \vec{u}
    // only initialization!
    rhoUFlux_
    (
        IOobject
        (
            "rhoUFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux_*linearInterpolate(U_)
    ),
    // fluxes in total energy equation: \varrho E \vec{u}
    // only initialization! 
    rhoEFlux_
    (
        IOobject
        (
            "rhoEFlux",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        rhoFlux_*linearInterpolate(he_ + 0.5*magSqr(U_) - (p_/rho_))
    ),
    dotX_
    (
        IOobject
        (
            "dotX",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimVelocity, vector::zero)
    ),
    gradp_
    (
        IOobject
        (
            "grad(p)",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", p_.dimensions()/dimLength, vector::zero)
    ),
    gradU_
    (
        IOobject
        (
            "grad(U)",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedTensor("zero", U_.dimensions()/dimLength, tensor::zero)
    ),
    gradrho_
    (
        IOobject
        (
            "grad(rho)",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", rho_.dimensions()/dimLength, vector::zero)
    ),
    gradk_
    (
        IOobject
        (
            "grad(k)",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector
        (
            "zero", 
            turbulenceModel.k()().dimensions()/dimLength, 
            vector::zero
        )
    ),
    coeffDict_(subOrEmptyDict(type + "Coeffs")),
    second_(lookup("secondOrder"))
{}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

autoPtr<godunovFlux> godunovFlux::New
(
    const volScalarField& rho,
    const volVectorField& U,
    const fluidThermo& thermoPhysicalModel,
    const turbulenceModel& turbulenceModel
)
{
    const word modelType
    (
        U.mesh().solutionDict().subDict("Riemann").lookup("fluxScheme")
    );

    Info<< "Selecting Godunov type flux scheme: " << modelType << endl;

    auto cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorIn
        (
            "godunovFlux::New"
            "("
                "const volVectorField&, "
                "const fluidThermo&, "
                "turbulenceModel&, "
            ")"
        )   << "Unknown godunovFlux type "
            << modelType << nl << nl
            << "Valid godunovFlux types:" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<godunovFlux>
    (
        cstrIter()(rho, U, thermoPhysicalModel, turbulenceModel)
    );
}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void godunovFlux::updateGrad()
{
    ITstream gScheme = lookup("reconGradScheme");

    gradp_ = fv::gradScheme<scalar>::New
        (
            mesh_, 
            gScheme
        )().grad(p_);

    const_cast<ITstream&>(gScheme).rewind();

    gradU_ = 
        fv::gradScheme<vector>::New
        (
            mesh_, 
            const_cast<ITstream&>(gScheme)
        )().grad(U_);

    const_cast<ITstream&>(gScheme).rewind();

    gradrho_ = 
        fv::gradScheme<scalar>::New
        (
            mesh_, 
            const_cast<ITstream&>(gScheme)
        )().grad(rho_);

    if (max(turbulenceModel_.k()).value() > 0.0)
    {
        const_cast<ITstream&>(gScheme).rewind();

        gradk_ = 
            fv::gradScheme<scalar>::New
            (
                mesh_, 
                const_cast<ITstream&>(gScheme)
            )().grad(turbulenceModel_.k());
    }
}


void godunovFlux::evaluateFlux()
{
   // Get face-to-cell addressing: face area point from owner to neighbour
    const labelUList& owner = mesh_.owner();
    const labelUList& neighbour = mesh_.neighbour();

    // Get the face area vector
    const surfaceVectorField& Sf = mesh_.Sf();
    const surfaceScalarField& magSf = mesh_.magSf();

    const volVectorField& C = mesh_.C();
    const surfaceVectorField& Cf = mesh_.Cf();

    const tmp<volScalarField> tk = turbulenceModel_.k();
    const volScalarField& k = tk();

    const volScalarField gamma
    (
        "gamma", 
        thermoPhysicalModel_.Cp()/thermoPhysicalModel_.Cv()
    );

    if (second_)
    {
        updateGrad();
    }

    // Calculate fluxes at internal faces
    forAll(owner, faceI)
    {
        label own = owner[faceI];
        label nei = neighbour[faceI];

        vector dL = Cf[faceI] - C[own];
        vector dR = Cf[faceI] - C[nei];

        // calculate fluxes with reconstructed primitive variables at faces
        // TODO: thermophysical variables are not reconstructed at faces!!!
        this->evaluate
        (
            rhoFlux_[faceI],
            rhoUFlux_[faceI],
            rhoEFlux_[faceI],
            Foam::max(p_[own] + second_*(dL & gradp_[own]), scalar(0)),
            Foam::max(p_[nei] + second_*(dR & gradp_[nei]), scalar(0)),
            U_[own] + second_*(dL & gradU_[own]),
            U_[nei] + second_*(dR & gradU_[nei]),
            rho_[own] + second_*(dL & gradrho_[own]),
            rho_[nei] + second_*(dR & gradrho_[nei]),
            k[own] + second_*(dL & gradk_[own]),
            k[nei] + second_*(dR & gradk_[nei]),
            gamma[own],       // left gamma
            gamma[nei],       // right gamma
            Sf[faceI],      // face vector
            magSf[faceI],   // face area
            dotX_[faceI]    // face velocity
        );
    }

   // Update boundary field and values
    forAll(p_.boundaryField(), patchI)
    {
        const labelUList& pOwner = mesh_.boundary()[patchI].faceCells();

        fvsPatchScalarField& pRhoFlux = rhoFlux_.boundaryFieldRef()[patchI];
        fvsPatchVectorField& pRhoUFlux = rhoUFlux_.boundaryFieldRef()[patchI];
        fvsPatchScalarField& pRhoEFlux = rhoEFlux_.boundaryFieldRef()[patchI];

        const fvPatchScalarField& pp = p_.boundaryField()[patchI];
        const fvPatchVectorField& pU = U_.boundaryField()[patchI];
        const fvPatchScalarField& prho = rho_.boundaryField()[patchI];
        const fvPatchScalarField& pk = k.boundaryField()[patchI];

        const fvPatchVectorField& pGradp = gradp_.boundaryField()[patchI];
        const fvPatchTensorField& pGradU = gradU_.boundaryField()[patchI];
        const fvPatchVectorField& pGradrho = gradrho_.boundaryField()[patchI];
        const fvPatchVectorField& pGradk = gradk_.boundaryField()[patchI];

        const fvPatchScalarField& pgamma = gamma.boundaryField()[patchI];

        const fvsPatchVectorField& pSf = Sf.boundaryField()[patchI];
        const fvsPatchScalarField& pMagSf = magSf.boundaryField()[patchI];
        const fvsPatchVectorField& pDotX = dotX_.boundaryField()[patchI];

        // special treatment at coupled boundaries
        if (pp.coupled())
        {
            // primitive variables
            const scalarField ppLeft(pp.patchInternalField());
            const scalarField ppRight(pp.patchNeighbourField());

            const vectorField pULeft(pU.patchInternalField());
            const vectorField pURight(pU.patchNeighbourField());

            const scalarField prhoLeft(prho.patchInternalField());
            const scalarField prhoRight(prho.patchNeighbourField());

            const scalarField pkLeft(pk.patchInternalField());
            const scalarField pkRight(pk.patchNeighbourField());

            const scalarField pgammaLeft(pgamma.patchInternalField());
            const scalarField pgammaRight(pgamma.patchNeighbourField());

            // cell gradients
            const vectorField pGradpLeft(pGradp.patchInternalField());
            const vectorField pGradpRight(pGradp.patchNeighbourField());

            const tensorField pGradULeft(pGradU.patchInternalField());
            const tensorField pGradURight(pGradU.patchNeighbourField());

            const vectorField pGradrhoLeft(pGradrho.patchInternalField());
            const vectorField pGradrhoRight(pGradrho.patchNeighbourField());

            const vectorField pGradkLeft(pGradk.patchInternalField());
            const vectorField pGradkRight(pGradk.patchNeighbourField());

            // cell and face centers
            vectorField pd(mesh_.Cf().boundaryField()[patchI].patch().delta());
            const vectorField& pCf = mesh_.Cf().boundaryField()[patchI];

            forAll(pOwner, faceI)
            {
                label own = pOwner[faceI];

                vector pdL = pCf[faceI] - C[own];
                vector pdR = pCf[faceI] - pd[faceI] - C[own];

                // Calculate fluxes at coupled boundary faces
                this->evaluate
                (
                    pRhoFlux[faceI],
                    pRhoUFlux[faceI],
                    pRhoEFlux[faceI],
                    Foam::max
                    (
                        ppLeft[faceI] + second_*(pdL & pGradpLeft[faceI]),
                        scalar(0)
                    ),
                    Foam::max
                    (
                        ppRight[faceI] + second_*(pdR & pGradpRight[faceI]),
                        scalar(0)
                    ),
                    pULeft[faceI] + second_*(pdL & pGradULeft[faceI]),
                    pURight[faceI] + second_*(pdR & pGradURight[faceI]),
                    prhoLeft[faceI] + second_*(pdL & pGradrhoLeft[faceI]),
                    prhoRight[faceI] + second_*(pdR & pGradrhoRight[faceI]),
                    pkLeft[faceI] + second_*(pdL & pGradkLeft[faceI]),
                    pkRight[faceI] + second_*(pdR & pGradkRight[faceI]),
                    pgammaLeft[faceI],  // face gamma
                    pgammaRight[faceI], // face gamma
                    pSf[faceI],         // face vector
                    pMagSf[faceI],      // face area
                    pDotX[faceI]        // face velocity
                );
            }
        }
        else
        {
            forAll(pp, faceI)
            {
                // Calculate fluxes at boundary faces
                this->evaluate
                (
                    pRhoFlux[faceI],
                    pRhoUFlux[faceI],
                    pRhoEFlux[faceI],
                    Foam::max(pp[faceI], scalar(0)),        // face p
                    Foam::max(pp[faceI], scalar(0)),        // face p
                    pU[faceI],        // face U
                    pU[faceI],        // face U
                    prho[faceI],      // face rho
                    prho[faceI],      // face rho
                    pk[faceI],        // face k
                    pk[faceI],        // face k
                    pgamma[faceI],    // face gamma
                    pgamma[faceI],    // face gamma
                    pSf[faceI],       // face vector
                    pMagSf[faceI],    // face area
                    pDotX[faceI]     // face velocity
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
