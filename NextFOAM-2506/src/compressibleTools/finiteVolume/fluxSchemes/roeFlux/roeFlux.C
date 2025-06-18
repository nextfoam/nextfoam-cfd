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

#include "roeFlux.H"
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

defineTypeNameAndDebug(roeFlux, 0);
addToRunTimeSelectionTable(godunovFlux, roeFlux, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

roeFlux::roeFlux
(
    const volScalarField& rho,
    const volVectorField& U,
    const fluidThermo& thermophysicalModel,
    const turbulenceModel& turbulenceModel,
    const word& modelName
)
:
    godunovFlux(modelName, rho, U, thermophysicalModel, turbulenceModel),
    epsilon_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "epsilon",
            coeffDict_,
            0.05
        )
    )
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void roeFlux::evaluate
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
    // Decode left and right total energy:
    const scalar E_L =
        pLeft/(rhoLeft*(gammaLeft -1.0)) + 0.5*magSqr(ULeft);
//      pLeft/(rhoLeft*(gammaLeft -1.0)) + 0.5*magSqr(ULeft) + kLeft; //cwo
    const scalar E_R =
        pRight/(rhoRight*(gammaRight-1.0)) + 0.5*magSqr(URight);
//      pRight/(rhoRight*(gammaRight-1.0)) + 0.5*magSqr(URight) + kRight; //cwo

    // normal vector
    vector normalVector = Sf/magSf;

    // Compute left and right contravariant velocities:
    const scalar contrVLeft = (ULeft & normalVector);
    const scalar contrVRight = (URight & normalVector);

    // Compute left and right total enthalpies:
    const scalar H_L = E_L + pLeft/rhoLeft;
    const scalar H_R = E_R + pRight/rhoRight;

    // Step 2: compute Roe averged quantities for face:
    const scalar rhoTilde = sqrt(max(rhoLeft*rhoRight,0.0));

    // Some temporary variables:
    const scalar rhoLeftSqrt = sqrt(max(rhoLeft,0.0));
    const scalar rhoRightSqrt = sqrt(max(rhoRight,0.0));

    const scalar wLeft = rhoLeftSqrt/(rhoLeftSqrt + rhoRightSqrt);
    const scalar wRight = 1.0 - wLeft;

    const vector UTilde = ULeft*wLeft + URight*wRight;
    const scalar hTilde = H_L*wLeft + H_R*wRight;
    const scalar qTildeSquare = magSqr(UTilde);
    const scalar gammaTilde = gammaLeft*wLeft + gammaRight*wRight;

    // Speed of sound
    const scalar cTilde =
        sqrt(max((gammaTilde - 1.0)*(hTilde - 0.5*qTildeSquare), 0.0));

    // Roe averaged contravariant velocity
    const scalar contrVTilde = (UTilde & normalVector) ;

    // Step 3: compute primitive differences:
    const scalar deltaP = pRight - pLeft;
    const scalar deltaRho = rhoRight - rhoLeft;
    const vector deltaU = URight - ULeft;
    const scalar deltaContrV = (deltaU & normalVector);

    // Step 4: compute wave strengths:

    // Roe and Pike - formulation
    const scalar r1 =
        (deltaP - rhoTilde*cTilde*deltaContrV)/(2.0*sqr(cTilde));
    const scalar r2 = deltaRho - deltaP/sqr(cTilde);
    const scalar r3 =
        (deltaP + rhoTilde*cTilde*deltaContrV)/(2.0*sqr(cTilde));

    // Step 5: compute right eigenvectors

    // rho row:
    const scalar l1rho = pTraits<scalar>::one;
    const scalar l2rho = pTraits<scalar>::one;
    const scalar l3rho = pTraits<scalar>::zero;
    const scalar l4rho = pTraits<scalar>::one;

    // first U column
    const vector l1U = UTilde - cTilde*normalVector;

    // second U column
    const vector l2U = UTilde;

    // third U column
    const vector l3U = deltaU - deltaContrV*normalVector;

    // fourth U column
    const vector l4U = UTilde + cTilde*normalVector;

    // E row
    const scalar l1e = hTilde - cTilde*contrVTilde;
    const scalar l2e = 0.5*qTildeSquare;
    const scalar l3e = (UTilde & deltaU) - contrVTilde*deltaContrV;
    const scalar l4e = hTilde + cTilde*contrVTilde;

    // Step 6: compute eigenvalues
    // derived from algebra by hand, only for Euler equation usefull
    const scalar Urot = dotX & normalVector;
    scalar lambda1 = mag(contrVTilde - cTilde - Urot);
    scalar lambda2 = mag(contrVTilde - Urot);
    scalar lambda3 = mag(contrVTilde + cTilde - Urot);

    // Step 7: check for Harten entropy correction
    //adjustable parameter
    const scalar eps = epsilon_.value()*max(lambda1,lambda3);

    if(lambda1 < eps)
    {
        lambda1 = (sqr(lambda1) + sqr(eps))/(2.0*eps);
    }

    if(lambda2 < eps)
    {
        lambda2 = (sqr(lambda2) + sqr(eps))/(2.0*eps);
    }

    if(lambda3 < eps)
    {
        lambda3 = (sqr(lambda3) + sqr(eps))/(2.0*eps);
    }

    // Step 8: Compute flux differences

    // Components of deltaF1
    const scalar diffF11 = lambda1*r1*l1rho;
    const vector diffF124 = lambda1*r1*l1U;
    const scalar diffF15 = lambda1*r1*l1e;

    // Components of deltaF2
    const scalar diffF21 = lambda2*(r2*l2rho + rhoTilde*l3rho);
    const vector diffF224 = lambda2*(r2*l2U + rhoTilde*l3U);
    const scalar diffF25 = lambda2*(r2*l2e + rhoTilde*l3e);

    // Components of deltaF3
    const scalar diffF31 = lambda3*r3*l4rho;
    const vector diffF324 = lambda3*r3*l4U;
    const scalar diffF35 = lambda3*r3*l4e;

    // Step 9: compute left and right fluxes

    // Left flux 5-vector
    const scalar fluxLeft11 = rhoLeft*contrVLeft;
    const vector fluxLeft124 = ULeft*fluxLeft11 + normalVector*pLeft;
    const scalar fluxLeft15 = H_L*fluxLeft11;

    // Right flux 5-vector
    const scalar fluxRight11 = rhoRight*contrVRight;
    const vector fluxRight124 = URight*fluxRight11 + normalVector*pRight;
    const scalar fluxRight15 = H_R*fluxRight11;

    // Step 10: compute face flux 5-vector
    const scalar flux1 =
        (
            fluxLeft11 + fluxRight11 
          - (diffF11 + diffF21 + diffF31)
          - (rhoLeft + rhoRight)*Urot
        )*0.5;

    const vector flux24 =
        (
            fluxLeft124 + fluxRight124 
          - (diffF124 + diffF224 + diffF324)
          - (rhoLeft*ULeft + rhoRight*URight)*Urot
        )*0.5;

    const scalar flux5 =
        (
            fluxLeft15  + fluxRight15 
          - (diffF15 + diffF25 + diffF35)
          - (rhoLeft*E_L + rhoRight*E_R)*Urot
        )*0.5;

    // Compute private data
    rhoFlux = flux1*magSf;
    rhoUFlux = flux24*magSf;
    rhoEFlux = flux5*magSf;
}


void roeFlux::update()
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
