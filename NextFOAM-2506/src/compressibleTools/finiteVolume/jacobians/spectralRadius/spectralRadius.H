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

Class
    Foam::compressible::spectralRadius

Description
    Viscous Jacobians are approximated by their spectral radii and linearized.

    References:
    Blazek, Jiri. "Computational Fluid Dynamics: Principles and Applications",
    2001, Referex Engineering, Elsevier, Chapter 6.2.4

Author
    JaeHeung Gill, NEXTfoam Co.,Ltd.

SourceFiles
    spectralRadius.C

\*---------------------------------------------------------------------------*/

#ifndef spectralRadius_H
#define spectralRadius_H

#include "fluidThermo.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace compressible
{

/*---------------------------------------------------------------------------*\
                           Class spectralRadius Declaration
\*---------------------------------------------------------------------------*/

class spectralRadius
{
    // Private data

        const Time& runTime_;
        const fvMesh& mesh_;

        const volScalarField& rho_;
        const volVectorField& U_;

        const volScalarField& p_;

        const fluidThermo& thermoPhysicalModel_;
        const turbulenceModel& turbulenceModel_;

        //- overrelaxation factor
        scalar w_;

        //- spectral radius of the inviscid flux Jacobian
        surfaceScalarField spectralRadiiInv_;

        //- spectral radius of the viscous flux Jacobian
        surfaceScalarField spectralRadiiVis_;

        //- spectral radius of the flux Jacobian from owner cells
        surfaceScalarField spectralRadiiown_;

        //- spectral radius of the flux Jacobian from neighbor cells
        surfaceScalarField spectralRadiinei_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        spectralRadius(const spectralRadius&);

        //- Disallow default bitwise assignment
        void operator=(const spectralRadius&);


public:

    // Constructors

        //- Construct from components
        spectralRadius
        (
            const volScalarField& rho,
            const volVectorField& U,
            const fluidThermo& thermoPhysicalModel,
            const turbulenceModel& turbulenceModel
        ); 


    // Destructor - default

    // Member Functions

        //- Return inviscid spectral radius
        const surfaceScalarField& spectralRadiiInv() const
        {
            return spectralRadiiInv_;
        }

        //- Return viscous spectral radius
        const surfaceScalarField& spectralRadiiVis() const
        {
            return spectralRadiiVis_;
        }

        //- Return spectral radius
        const surfaceScalarField& spectralRadiiown() const
        {
            return spectralRadiiown_;
        }

        //- Return spectral radius
        const surfaceScalarField& spectralRadiinei() const
        {
            return spectralRadiinei_;
        }

        // Evaluate spectral radius
        void evaluateSpectralRadii
        (
            scalar& spectralRadiiInv,
            scalar& spectralRadiiVis,
            scalar& spectralRadiiown,
            scalar& spectralRadiinei,

            vector& dLoc,
            const scalar alphaEffLeft,
            const scalar alphaEffRight,

            const scalar pLeft,
            const scalar pRight,
            const vector ULeft,
            const vector URight,
            const scalar rhoLeft,
            const scalar rhoRight,

            const scalar gammaLeft,
            const scalar gammaRight,

            const vector Sf,
            const scalar magSf
        ) const; 

        //- Update spectral radius 
        void update();

        //- Return sum of inviscid eigenvalue for each cell
        const tmp<scalarField> LambdaC() const;

        //- Return sum of viscous eigenvalue for each cell
        const tmp<scalarField> LambdaV() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
