/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2017 OpenCFD Ltd.
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
    Foam::multiRegionSimpleControl

Description
    SIMPLE control class to supply convergence information/checks for
    the multi-region SIMPLE loop.

\*---------------------------------------------------------------------------*/

#ifndef multiRegionSimpleControl_H
#define multiRegionSimpleControl_H

#include "multiRegionSolutionControl.H"

//- Declare that SIMPLE Control will be used
#define SIMPLE_CONTROL
#define CREATE_MESH createMultiRegionMeshPostProcess.H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class multiRegionSimpleControl Declaration
\*---------------------------------------------------------------------------*/

class multiRegionSimpleControl
:
    public multiRegionSolutionControl
{

protected:

    // Protected data


        // Solution controls

            //- Initialised flag
            bool initialised_;


    // Protected Member Functions

        //- Read controls from fvSolution dictionary
        void readControl();

        //- Return true if all convergence checks are satisfied
        bool criteriaSatisfied();


private:

    // Private member functions

        //- No copy construct
        multiRegionSimpleControl(const multiRegionSimpleControl&) = delete;

        //- No copy assignment
        void operator=(const multiRegionSimpleControl&) = delete;


public:

    // Static Data Members

        //- Run-time type information
        TypeName("multiRegionSimpleControl");


    // Constructors

        //- Construct from mesh and the name of control sub-dictionary
        multiRegionSimpleControl
        (
            const Time& runTime, 
            const PtrList<fvMesh>& fluid,
            const PtrList<fvMesh>& solid,
            bool foamLog = true
        );


    //- Destructor
    virtual ~multiRegionSimpleControl() = default;


    // Member Functions

        // Solution control

            //- SIMPLE loop
            bool loop();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
