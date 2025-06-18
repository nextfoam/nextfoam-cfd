/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "hllcFlux.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace compressible
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace godunovFluxes
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(hllcFlux, 0);
addToRunTimeSelectionTable(godunovFlux, hllcFlux, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

hllcFlux::hllcFlux
(
    const volScalarField& rho,
    const volVectorField& U,
    const fluidThermo& thermophysicalModel,
    const turbulenceModel& turbulenceModel,
    const word& modelName
)
:
    godunovFlux(modelName, rho, U, thermophysicalModel, turbulenceModel)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void hllcFlux::evaluate
(   
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    const scalar pLeft,
    const scalar pRight,
    const vector ULeft,
    const vector URight,
    const scalar rhoLeft,
    const scalar rhoRight,
    const scalar kLeft,
    const scalar kRight,
    const scalar gammaLeft,
    const scalar gammaRight,
    const vector Sf,
    const scalar magSf,
    const vector dotX
) const
{
    // bounding variables
    const scalar rhoMin = SMALL;

    // Step 1: decode left and right:
    // normal vector
    const vector normalVector = Sf/magSf;

    // speed of sound, for left and right side, assuming perfect gas
    const scalar aLeft = 
        Foam::sqrt(max(SMALL, gammaLeft*pLeft/max(rhoLeft, rhoMin)));
    const scalar aRight = 
        Foam::sqrt(max(SMALL, gammaRight*pRight/max(rhoRight, rhoMin)));
//    const scalar aLeft  = Foam::sqrt(gammaLeft *pLeft /rhoLeft);
//      const scalar aRight = Foam::sqrt(gammaRight*pRight/rhoRight);

    // DensityTotalEnergy
    const scalar rhoELeft =
        pLeft/(gammaLeft - 1.0) + rhoLeft*(0.5*magSqr(ULeft) + kLeft);
    const scalar rhoERight =
        pRight/(gammaRight - 1.0) + rhoRight*(0.5*magSqr(URight) + kRight);

    // DensityVelocity
    const vector rhoULeft  = rhoLeft*ULeft;
    const vector rhoURight = rhoRight*URight;

    // Compute left and right total enthalpies:
    const scalar HLeft  = (rhoELeft + pLeft)/max(rhoLeft, rhoMin);
    const scalar HRight = (rhoERight + pRight)/max(rhoRight, rhoMin);

    // compute velocity relative to mesh
    const vector URelLeft  = ULeft - dotX;
    const vector URelRight = URight - dotX;

    // compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
    const scalar qLeft = (URelLeft & normalVector);
    const scalar qRight = (URelRight & normalVector);

    // Step 2:
    // needs rho_{l,r}, U_{l,r}, H_{l,r}, gamma_{l,r}, Gamma_{l,r}, q_{l,r}
    // compute Roe weights
    const scalar rhoLeftSqrt = Foam::sqrt(max(rhoLeft, rhoMin));
    const scalar rhoRightSqrt = Foam::sqrt(max(rhoRight, rhoMin));

    const scalar wLeft = rhoLeftSqrt/max((rhoLeftSqrt + rhoRightSqrt), rhoMin);
    const scalar wRight = 1.0 - wLeft;

    // Roe averaged velocity
    const vector UTilde = wLeft*ULeft + wRight*URight;

    // Roe averaged relative face velocity
    const scalar contrURelTilde = ((UTilde - dotX) & normalVector);

    // Roe averaged contravariant velocity
    const scalar contrUTilde = (UTilde & normalVector);

    // Roe averaged total enthalpy
    const scalar HTilde = wLeft*HLeft + wRight*HRight;

    // Roe averaged gamma
    // TODO: needs to be verified for non constant gamma values
    const scalar gammaTilde = wLeft*gammaLeft + wRight*gammaRight;

    // Speed of sound with Roe reconstruction values
    const scalar aTilde
    (
        Foam::sqrt
        (
            max(SMALL, (gammaTilde - 1.0)*(HTilde - 0.5*sqr(contrUTilde)))
        )
    );

    // Step 3: compute signal speeds for face:
    const scalar SLeft = min(qLeft - aLeft,  contrURelTilde - aTilde);
    const scalar SRight = max(qRight + aRight, contrURelTilde + aTilde);

    // See for an alternative computation of wave speed: Toro (1999) 
    // Eq. 10.46-10.48

    const scalar SStar
    (
        (
            rhoRight*qRight*(SRight - qRight) 
          - rhoLeft*qLeft*(SLeft - qLeft) 
          + pLeft - pRight 
        )
        /
        stabilise
        (
            (rhoRight*(SRight - qRight) - rhoLeft*(SLeft - qLeft)),
            SMALL
        )
    );

    // compute pressure in star region from the right side
    const scalar pStarRight =
        rhoRight*(qRight - SRight)*(qRight - SStar) + pRight;

    // should be equal to the left side
    const scalar pStarLeft  =
        rhoLeft*(qLeft -  SLeft)*(qLeft - SStar) + pLeft;

    // give a warning if this is not the case
    if ( mag(pStarRight - pStarLeft) > 1e-6 )
    {
            Info << "mag(pStarRight - pStarLeft) > SMALL " << endl;
    }

    // use pStarRight for pStar, as in theory, pStarRight == pStarLeft
    const scalar pStar = pStarRight;

    // Step 4: upwinding - compute states:
    scalar convectionSpeed = 0.0;
    scalar rhoState = 0.0;
    vector rhoUState = vector::zero;
    scalar rhoEState = 0.0;
    scalar pState = 0.0;

    // TODO: Maybe one can use pos/neg implementation, but then one has to
    // evaluate all 4 states at each iteration!
    // label A = pos(SLeft);
    // label B = pos(SStar);
    // label C = pos(SRight);
    // please double check the bool operators again, if one want's to
    // implement this!!!
    // scalar convectionSpeed = A*B*C*qLeft+(1-A)*B*C*SStar
    //     +(1-A)*(1-B)*C*SStar+(1-A)*(1-B)*(1-C)*qRight:

    if (pos(SLeft))
    {
        // compute F_l
        convectionSpeed = qLeft;
        rhoState = rhoLeft;
        rhoUState = rhoULeft;
        rhoEState = rhoELeft;
        pState = pLeft;
    }
    else if (pos(SStar))
    {
        scalar omegaLeft = scalar(1.0)/stabilise((SLeft - SStar),SMALL);

        // compute left star region
        convectionSpeed = SStar;
        rhoState = omegaLeft*(SLeft - qLeft)*rhoLeft;
        rhoUState =  
        (
            (
                (SLeft - qLeft)*rhoULeft + (pStar - pLeft)*normalVector
            )*omegaLeft
        );
        rhoEState = 
        (
            (
                (SLeft - qLeft)*rhoELeft - pLeft*qLeft + pStar*SStar
            )*omegaLeft
        );
        pState = pStar;
    }
    else if (pos(SRight))
    {
        scalar omegaRight = scalar(1.0)/stabilise((SRight - SStar), SMALL);

        // compute right star region
        convectionSpeed = SStar;
        rhoState = omegaRight*(SRight - qRight)*rhoRight;
        rhoUState =
        (  
            (
                (SRight - qRight)*rhoURight + (pStar - pRight)*normalVector
            )*omegaRight
        );
        rhoEState = 
        (
            (
                (SRight - qRight)*rhoERight - pRight*qRight + pStar*SStar 
            )*omegaRight
        );
        pState = pStar;
    }
    else if (neg(SRight))
    {
        // compute F_r
        convectionSpeed = qRight;
        rhoState = rhoRight;
        rhoUState = rhoURight;
        rhoEState = rhoERight;
        pState = pRight;
    }
    else
    {
        Info << "Error in HLLC Riemann solver" << endl;
    }

    rhoFlux  = (convectionSpeed*rhoState)*magSf;
    rhoUFlux = (convectionSpeed*rhoUState + pState*normalVector)*magSf;
    rhoEFlux = 
    (
        (
            convectionSpeed*(rhoEState + pState) 
          + pState*(dotX & normalVector)
        )*magSf
    );
}


void hllcFlux::update()
{
	evaluateFlux();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace godunovFluxes

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
