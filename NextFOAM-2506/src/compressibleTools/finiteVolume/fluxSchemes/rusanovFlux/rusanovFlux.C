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

#include "rusanovFlux.H"
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

defineTypeNameAndDebug(rusanovFlux, 0);
addToRunTimeSelectionTable(godunovFlux, rusanovFlux, dictionary);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rusanovFlux::rusanovFlux
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

void rusanovFlux::evaluate
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
//    const scalar aRight = Foam::sqrt(gammaRight*pRight/rhoRight);

    // DensityTotalEnergy
    const scalar rhoELeft =
        pLeft/(gammaLeft - 1.0) + rhoLeft*(0.5*magSqr(ULeft) + kLeft);
    const scalar rhoERight =
        pRight/(gammaRight - 1.0) + rhoRight*(0.5*magSqr(URight) + kRight);

    // DensityVelocity
    const vector rhoULeft  = rhoLeft*ULeft;
    const vector rhoURight = rhoRight*URight;

    // Compute left and right total enthalpies:
//     const scalar rhoHLeft  = (rhoELeft  + pLeft);
//     const scalar rhoHRight = (rhoERight + pRight);

    // compute qLeft and qRight (q_{l,r} = U_{l,r} \bullet n)
    const scalar qLeft = ((ULeft  - dotX) & normalVector);
    const scalar qRight = ((URight - dotX) & normalVector);

    const scalar lambda1 = mag(qLeft) + aLeft;
    const scalar lambda2 = mag(qRight) + aRight;

//     const scalar lambda1 = mag(ULeft  & normalVector) + aLeft;
//     const scalar lambda2 = mag(URight & normalVector) + aRight;

    const scalar lambdaMax = max(lambda1, lambda2);

    // Step 10: compute face flux 5-vector
    rhoFlux =
    (
        (
            qLeft*rhoLeft + qRight*rhoRight - lambdaMax*(rhoRight-rhoLeft)
        )*magSf/2.0
    );

    rhoUFlux =
    (
        (
            qLeft*rhoULeft + qRight*rhoURight 
          + (pLeft + pRight)*normalVector 
          - lambdaMax*(rhoURight - rhoULeft)
        )*magSf/2.0
    );

    rhoEFlux =
    (
        (
            qLeft*rhoELeft + qRight*rhoERight 
          + ((pLeft*ULeft + pRight*URight) & normalVector) 
          - lambdaMax*(rhoERight - rhoELeft)
        )*magSf/2.0
    );


// Backup
/*

    // Adiabatic exponent is constant for ideal gas but if Cp=Cp(T)
    // it must be computed for each cell and evaluated at each face
    // through reconstruction
//     const scalar kappaLeft = (CvLeft+RLeft)/CvLeft;
//     const scalar kappaRight = (CvRight+RRight)/CvRight;

    // Step 2: compute Roe averged quantities for face:
//     const scalar rhoTilde = sqrt(max(rhoLeft*rhoRight,0.0));

    // Some temporary variables:
//     const scalar rhoLeftSqrt = sqrt(max(rhoLeft,0.0));
//     const scalar rhoRightSqrt = sqrt(max(rhoRight,0.0));

//     const scalar wLeft = rhoLeftSqrt/(rhoLeftSqrt + rhoRightSqrt);
//     const scalar wRight = 1.0 - wLeft;

//     const vector UTilde = ULeft*wLeft + URight*wRight;
//     const scalar hTilde = hLeft*wLeft + hRight*wRight;
//     const scalar qTildeSquare = magSqr(UTilde);
//     const scalar kappaTilde = kappaLeft*wLeft + kappaRight*wRight;

//     // Speed of sound
//     const scalar cTilde =
//         sqrt(max((kappaTilde - 1.0)*(hTilde - 0.5*qTildeSquare),0.0));

    // Roe averaged contravariant velocity
//     const scalar contrVTilde = (UTilde & normalVector) ;

    // Step 3: compute primitive differences:
//     const scalar deltaP = pRight - pLeft;
//     const scalar deltaRho = rhoRight - rhoLeft;
//     const vector deltaU = URight - ULeft;
//     const scalar deltaContrV = (deltaU & normalVector);

    // Step 4: compute wave strengths:

//     // Roe and Pike - formulation
//     const scalar r1 =
//         (deltaP - rhoTilde*cTilde*deltaContrV)/(2.0*sqr(cTilde));
//     const scalar r2 = deltaRho - deltaP/sqr(cTilde);
//     const scalar r3 =
//         (deltaP + rhoTilde*cTilde*deltaContrV)/(2.0*sqr(cTilde));

    // Step 5: compute l vectors

//     // rho row:
//     const scalar l1rho = pTraits<scalar>::one;
//     const scalar l2rho = pTraits<scalar>::one;
//     const scalar l3rho = pTraits<scalar>::zero;
//     const scalar l4rho = pTraits<scalar>::one;
// 
//     // first U column
//     const vector l1U = UTilde - cTilde*normalVector;
// 
//     // second U column
//     const vector l2U = UTilde;
// 
//     // third U column
//     const vector l3U = deltaU - deltaContrV*normalVector;
// 
//     // fourth U column
//     const vector l4U = UTilde + cTilde*normalVector;
// 
//     // E row
//     const scalar l1e = hTilde - cTilde*contrVTilde;
//     const scalar l2e = 0.5*qTildeSquare;
//     const scalar l3e = (UTilde & deltaU) - contrVTilde*deltaContrV;
//     const scalar l4e = hTilde + cTilde*contrVTilde;

    // Step 6: compute eigenvalues
    // derived from algebra by hand, only for Euler equation usefull

    const scalar Urot = dotX & normalVector;
//     scalar lambda1 = mag(contrVTilde - cTilde - Urot);
//     scalar lambda2 = mag(contrVTilde - Urot);
//     scalar lambda3 = mag(contrVTilde + cTilde - Urot);
// 
//     scalar lambdaMax = max(max(lambda1,lambda2),lambda3);
//     scalar lambdaMax = mag(contrVTilde - Urot)+cTilde;

    // speed of sound, for left and right side, assuming perfect gas
    const scalar aLeft  = Foam::sqrt(max(SMALL,(RLeft /CvLeft +1.0)*RLeft *TLeft));
    const scalar aRight = Foam::sqrt(max(SMALL,(RRight/CvRight+1.0)*RRight*TRight));

    const scalar lambda1 = mag((ULeft -dotX) & normalVector) + aLeft;
    const scalar lambda2 = mag((URight-dotX) & normalVector) + aRight;

//     const scalar lambda1 = mag((ULeft)  & normalVector) + aLeft;
//     const scalar lambda2 = mag((URight) & normalVector) + aRight;

    const scalar lambdaMax = max(lambda1,lambda2);

    // Step 7: Compute flux differences

//     // Components of deltaF1
//     const scalar diffF11  = lambdaMax*r1*l1rho;
//     const vector diffF124 = lambdaMax*r1*l1U;
//     const scalar diffF15  = lambdaMax*r1*l1e;
// 
//     // Components of deltaF2
//     const scalar diffF21  = lambdaMax*(r2*l2rho + rhoTilde*l3rho);
//     const vector diffF224 = lambdaMax*(r2*l2U   + rhoTilde*l3U);
//     const scalar diffF25  = lambdaMax*(r2*l2e   + rhoTilde*l3e);
// 
//     // Components of deltaF3
//     const scalar diffF31  = lambdaMax*r3*l4rho;
//     const vector diffF324 = lambdaMax*r3*l4U;
//     const scalar diffF35  = lambdaMax*r3*l4e;

    // Step 8: compute left and right fluxes

    // Left flux 5-vector
    const scalar fluxLeft11  = rhoLeft*contrVLeft;
    const vector fluxLeft124 = ULeft*fluxLeft11 + normalVector*pLeft;
    const scalar fluxLeft15  = hLeft*fluxLeft11;

    // Right flux 5-vector
    const scalar fluxRight11  = rhoRight*contrVRight;
    const vector fluxRight124 = URight*fluxRight11 + normalVector*pRight;
    const scalar fluxRight15  = hRight*fluxRight11;

//     // Step 10: compute face flux 5-vector
//     const scalar flux1 =
//         0.5*(fluxLeft11  + fluxRight11  -(rhoLeft+rhoRight)*Urot              - (diffF11  + diffF21  + diffF31));
//     const vector flux24 =
//         0.5*(fluxLeft124 + fluxRight124 -(rhoLeft*ULeft+rhoRight*URight)*Urot - (diffF124 + diffF224 + diffF324));
//     const scalar flux5 =
//         0.5*(fluxLeft15  + fluxRight15  -(rhoLeft*eLeft+rhoRight*eRight)*Urot - (diffF15  + diffF25  + diffF35));

    // Step 10: compute face flux 5-vector
    const scalar flux1 =
        0.5*(fluxLeft11  + fluxRight11  -(rhoLeft+rhoRight)*Urot              - lambdaMax*(rhoRight-rhoLeft));
    const vector flux24 =
        0.5*(fluxLeft124 + fluxRight124 -(rhoLeft*ULeft+rhoRight*URight)*Urot - lambdaMax*(rhoRight*URight-rhoLeft*ULeft));
    const scalar flux5 =
        0.5*(fluxLeft15  + fluxRight15  -(rhoLeft*eLeft+rhoRight*eRight)*Urot - lambdaMax*(rhoRight*eRight-rhoLeft*eLeft));

    // Compute private data
    rhoFlux  = flux1 *magSf;
    rhoUFlux = flux24*magSf;
    rhoEFlux = flux5 *magSf;
*/
}


void rusanovFlux::update()
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
