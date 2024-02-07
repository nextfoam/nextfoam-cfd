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

#include "AUSMplusFlux.H"
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

defineTypeNameAndDebug(AUSMplusFlux, 0);
addToRunTimeSelectionTable(godunovFlux, AUSMplusFlux, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

AUSMplusFlux::AUSMplusFlux
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

void AUSMplusFlux::evaluate
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
    const scalar alpha = 3.0/16.0;
    const scalar beta = 1.0/8.0;

    // bounding variables
    const scalar rhoMin = SMALL;

    // normal vector
    const vector normalVector = Sf/magSf;

    // speed of sound, for left and right side, assuming perfect gas
    const scalar aLeft =
        Foam::sqrt(max(SMALL, gammaLeft*pLeft/max(rhoLeft, rhoMin)));
    const scalar aRight =
        Foam::sqrt(max(SMALL, gammaRight*pRight/max(rhoRight, rhoMin)));

    // DensityTotalEnergy
    const scalar rhoHLeft
    (
        pLeft*gammaLeft/(gammaLeft - 1.0) 
      + rhoLeft*(0.5*magSqr(ULeft) + kLeft)
    );
    const scalar rhoHRight
    (
        pRight*gammaRight/(gammaRight - 1.0) 
      + rhoRight*(0.5*magSqr(URight) + kRight)
    );

    // Step 1: decode left and right:

    // Compute conservative variables assuming perfect gas law

    // DensityVelocity
    const vector rhoULeft = rhoLeft *ULeft;
    const vector rhoURight = rhoRight*URight;

    // compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
    const scalar qLeft = ((ULeft - dotX) & normalVector);
    const scalar qRight = ((URight - dotX) & normalVector);

    const scalar aTilde = 0.5*(aLeft + aRight);

    const scalar MaRelLeft = qLeft/aTilde;
    const scalar MaRelRight = qRight/aTilde;

    const scalar magMaRelLeft = mag(MaRelLeft);
    const scalar magMaRelRight = mag(MaRelRight);

    const scalar Ma1PlusLeft = 0.5*(MaRelLeft + magMaRelLeft);
    const scalar Ma1MinusRight = 0.5*(MaRelRight - magMaRelRight);

    const scalar Ma2PlusLeft =  0.25*sqr(MaRelLeft + 1.0);
    const scalar Ma2PlusRight =  0.25*sqr(MaRelRight + 1.0);
    const scalar Ma2MinusLeft = -0.25*sqr(MaRelLeft - 1.0);
    const scalar Ma2MinusRight = -0.25*sqr(MaRelRight - 1.0);

    const scalar Ma4BetaPlusLeft
    (
        (magMaRelLeft >= 1.0) 
      ? Ma1PlusLeft   
      : Ma2PlusLeft*(1.0 - 16.0*beta*Ma2MinusLeft)
    );

    const scalar Ma4BetaMinusRight
    (
        (magMaRelRight >= 1.0) 
      ? Ma1MinusRight 
      : Ma2MinusRight*(1.0 + 16.0*beta*Ma2PlusRight)
    );

    const scalar P5alphaPlusLeft
    (
        (magMaRelLeft >= 1.0) 
      ? Ma1PlusLeft/MaRelLeft    
      : Ma2PlusLeft*(2.0 - MaRelLeft - 16.0*alpha*MaRelLeft*Ma2MinusLeft)
    );

    const scalar P5alphaMinusRight
    (
        (magMaRelRight >= 1.0) 
      ? Ma1MinusRight/MaRelRight 
      : Ma2MinusRight*(-2.0 - MaRelRight + 16.0*alpha*MaRelRight*Ma2PlusRight)
    );

    const scalar MaRelTilde = Ma4BetaPlusLeft + Ma4BetaMinusRight;
    const scalar pTilde = pLeft*P5alphaPlusLeft + pRight*P5alphaMinusRight;

    const scalar URelTilde = MaRelTilde*aTilde;
    const scalar magURelTilde = mag(MaRelTilde)*aTilde;
    // There is a typo in Luo et. al, J. Comp. Physics 194 (2004), 
    // Chap 4.2 Eq. 4.8
    // refer to the origial Paper from Liou, J. Comp. Physics 129 (1996), 
    // Chap4, Eq. 42
    rhoFlux = 
    (
        (
            URelTilde*(rhoLeft + rhoRight) 
          - magURelTilde*(rhoRight - rhoLeft)
        )*magSf/2.0
    );

    rhoUFlux = 
    (
        (
            (
                URelTilde*(rhoULeft + rhoURight)
              - magURelTilde*(rhoURight - rhoULeft)
            )/2.0 + pTilde*normalVector
        )*magSf
    );

    rhoEFlux = 
    (
        (
            (
                URelTilde*(rhoHLeft + rhoHRight)
              - magURelTilde*(rhoHRight - rhoHLeft)
            )/2.0 + pTilde*(dotX & normalVector)
        )*magSf
    );
}


void AUSMplusFlux::update()
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
