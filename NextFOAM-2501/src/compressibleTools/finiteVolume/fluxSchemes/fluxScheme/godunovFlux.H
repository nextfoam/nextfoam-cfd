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

Class
    Foam::compressible::godunovFlux

Description
    Generic Godunov flux class. Does an extrapolation of cell centered
    variables to face centeres, using a Taylor-Series Expansion. Therefore a
    slope limiter is utilized.

    References:
    Weiss, Jonathan M. and Wayne A. Smith, "Preconditioning Applied to Variable
    and Constant Density Flows." AIAA Journal 33/11 (November 1995): 2050–2057.

Author
    Oliver Borm  All rights reserved.
    Modifiled JaeHeung Gill, NEXTFoam Co.,Ltd.

SourceFiles
    godunovFlux.C

\*---------------------------------------------------------------------------*/

#ifndef godunovFlux_H
#define godunovFlux_H

#include "fluidThermo.H"
#include "turbulentFluidThermoModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace compressible
{

/*---------------------------------------------------------------------------*\
                           Class godunovFlux Declaration
\*---------------------------------------------------------------------------*/

class godunovFlux
:
    public dictionary
{

protected:

    // Protected data
    
        const Time& runTime_;
        const fvMesh& mesh_;

        const volScalarField& rho_;
        const volVectorField& U_;

        const volScalarField& p_;
        const volScalarField& he_;

        const fluidThermo& thermoPhysicalModel_;
        const turbulenceModel& turbulenceModel_;

        surfaceScalarField rhoFlux_;
        surfaceVectorField rhoUFlux_;
        surfaceScalarField rhoEFlux_;

        //- mesh fluxes
        surfaceVectorField dotX_;

        //- static Pressure
        volVectorField gradp_;

        //- Velocity
        volTensorField gradU_;

        //- Density
        volVectorField gradrho_;
 
        //- turbulent kinetic energy
        volVectorField gradk_;

        //- Model coefficients dictionary
        dictionary coeffDict_;

		Switch second_;


private:

    // Protected Member Functions

        //- Disallow default bitwise copy construct
        godunovFlux(const godunovFlux&);

        //- Disallow default bitwise assignment
        void operator=(const godunovFlux&);


public:

    //- Runtime type information
    TypeName("godunovFlux");


    // Declare run-time constructor selection table

        declareRunTimeSelectionTable
        (
            autoPtr,
            godunovFlux,
            dictionary,
            (
                const volScalarField& rho,
                const volVectorField& U,
                const fluidThermo& thermoPhysicalModel,
                const turbulenceModel& turbulenceModel
            ),
            (rho, U, thermoPhysicalModel, turbulenceModel)
        );


    // Constructors

        //- Construct from components
        godunovFlux
        (
            const word& type,
            const volScalarField& rho,
            const volVectorField& U,
            const fluidThermo& thermoPhysicalModel,
            const turbulenceModel& turbulenceModel
        );


    // Selectors

        //- Return a reference to the selected godunovFlux model
        static autoPtr<godunovFlux> New
        (
            const volScalarField& rho,
            const volVectorField& U,
            const fluidThermo& thermoPhysicalModel,
            const turbulenceModel& turbulenceModel
        );


    //- Destructor
    virtual ~godunovFlux()
    {}


    // Member Functions

        //- Return density flux
        const surfaceScalarField& rhoFlux() const
        {
            return rhoFlux_;
        }

        //- Return velocity flux
        const surfaceVectorField& rhoUFlux() const
        {
            return rhoUFlux_;
        }

        //- Return energy flux
        const surfaceScalarField& rhoEFlux() const
        {
            return rhoEFlux_;
        }

        //- Return mesh flux
        const surfaceVectorField& dotX() const
        {
            return dotX_;
        }

        //- Return mesh flux
        surfaceVectorField& dotX()
        {
            return dotX_;
        }

        void updateGrad();

        void evaluateFlux();

        virtual void evaluate
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
        ) const = 0;

        //- Update fluxes based on current state
        virtual void update() = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace compressible

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
