/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2018-2021 OpenCFD Ltd.
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

#include "signalHandler.H"
#include "Pstream.H"
#include "OSspecific.H"
#include "IFstream.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::fileName Foam::signalHandler::procInfoFile_ = "";


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::signalHandler::updateProcState()
{
    if (Pstream::master())
    {
        if (Foam::isFile(procInfoFile_))
        {
            IFstream is(procInfoFile_);
            dictionary procInfo(is);
            procInfo.add("State", "stopped", true);
            OFstream os(procInfoFile_);
            procInfo.writeEntries(os, false);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::signalHandler::signalHandler()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::signalHandler::~signalHandler()
{}


// ************************************************************************* //
