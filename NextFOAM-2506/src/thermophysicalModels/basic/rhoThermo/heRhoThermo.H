/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2017 OpenCFD Ltd.
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
    Foam::heRhoThermo

Description
    Energy for a mixture based on density

SourceFiles
    heRhoThermo.C

\*---------------------------------------------------------------------------*/

#ifndef heRhoThermo_H
#define heRhoThermo_H

#include "rhoThermo.H"
#include "heThermo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class heRhoThermo Declaration
\*---------------------------------------------------------------------------*/

template<class BasicPsiThermo, class MixtureType>
class heRhoThermo
:
    public heThermo<BasicPsiThermo, MixtureType>
{
    // Private Member Functions

        //- Calculate the thermo variables
        void calculate
        (
            const volScalarField& p,
            volScalarField& T,
            volScalarField& he,
            volScalarField& psi,
            volScalarField& rho,
            volScalarField& mu,
            volScalarField& alpha,
            const bool doOldTimes
        );

        //- Calculate the rho and mu only
        void calculateRhoMu
        (
            const volScalarField& p,
            volScalarField& T,
            volScalarField& rho,
            volScalarField& mu,
            const bool doOldTimes
        );

        //- Calculate the internal thermo variables
        void calculateInternal
        (
            const volScalarField& p,
            volScalarField& T,
            volScalarField& he,
            volScalarField& psi,
            volScalarField& rho,
            volScalarField& mu,
            volScalarField& alpha,
            const bool doOldTimes
        );

        //- Calculate the boundary thermo variables
        void calculateBoundary
        (
            const volScalarField& p,
            volScalarField& T,
            volScalarField& he,
            volScalarField& psi,
            volScalarField& rho,
            volScalarField& mu,
            volScalarField& alpha,
            const bool doOldTimes
        );

        void addHydrostaticPressure
        (
            volScalarField& p,
            const volScalarField& rho
        ); // by Gill

        //- Construct as copy (not implemented)
        heRhoThermo(const heRhoThermo<BasicPsiThermo, MixtureType>&);


public:

    //- Runtime type information
    TypeName("heRhoThermo");


    // Constructors

        //- Construct from mesh and phase name
        heRhoThermo
        (
            const fvMesh&,
            const word& phaseName
        );


        //- Construct from mesh,dictionary,phase name with a single temperature
        heRhoThermo
        (
            const fvMesh&,
            const word& phaseName,
            const word& dictName
        );


    //- Destructor
    virtual ~heRhoThermo();


    // Member functions

        //- Update properties
        virtual void correct();

        //- Update rho and mu only
        virtual void correctRhoMu();

        //- Update internal field properties
        virtual void correctInternal();

        //- Update boundary field properties
        virtual void correctBoundary();

        //- Update properties from T
        //virtual void correctT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "heRhoThermo.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
