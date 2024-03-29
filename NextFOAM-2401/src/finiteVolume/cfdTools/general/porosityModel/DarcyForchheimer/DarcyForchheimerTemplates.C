/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2012-2016 OpenFOAM Foundation
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

#include "volFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class RhoFieldType>
void Foam::porosityModels::DarcyForchheimer::apply
(
    scalarField& Udiag,
    vectorField& Usource,
    const scalarField& V,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{
    if (mesh_.foundObject<volVectorField>("porousGradp"))
    {
        volVectorField& porousGradp
        (
            const_cast<volVectorField&>
            (
                mesh_.lookupObject<volVectorField>("porousGradp")
            )
        );

        forAll(cellZoneIDs_, zoneI)
        {
            const tensorField& dZones = D_[zoneI];
            const tensorField& fZones = F_[zoneI];

            const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

            forAll(cells, i)
            {
                const label celli = cells[i];
                const label j = this->fieldIndex(i);
                const tensor Cd =
                    mu[celli]*dZones[j] + (rho[celli]*mag(U[celli]))*fZones[j];

                const scalar isoCd = tr(Cd);

                Udiag[celli] += V[celli]*isoCd;
                Usource[celli] -= V[celli]*((Cd - I*isoCd) & U[celli]);

                //porousGradp.primitiveFieldRef()[celli] = -(Cd &  U[celli]);
            }
        }

        //porousGradp.correctBoundaryConditions();
    }
    else
    {
        forAll(cellZoneIDs_, zoneI)
        {
            const tensorField& dZones = D_[zoneI];
            const tensorField& fZones = F_[zoneI];

            const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

            forAll(cells, i)
            {
                const label celli = cells[i];
                const label j = this->fieldIndex(i);
                const tensor Cd =
                    mu[celli]*dZones[j] + (rho[celli]*mag(U[celli]))*fZones[j];

                const scalar isoCd = tr(Cd);

                Udiag[celli] += V[celli]*isoCd;
                Usource[celli] -= V[celli]*((Cd - I*isoCd) & U[celli]);
            }
        }
    }
}


template<class RhoFieldType>
void Foam::porosityModels::DarcyForchheimer::apply
(
    tensorField& AU,
    const RhoFieldType& rho,
    const scalarField& mu,
    const vectorField& U
) const
{
    if (mesh_.foundObject<volVectorField>("porousGradp"))
    {
        volVectorField& porousGradp
        (
            const_cast<volVectorField&>
            (
                mesh_.lookupObject<volVectorField>("porousGradp")
            )
        );

        forAll(cellZoneIDs_, zoneI)
        {
            const tensorField& dZones = D_[zoneI];
            const tensorField& fZones = F_[zoneI];

            const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

            forAll(cells, i)
            {
                const label celli = cells[i];
                const label j = this->fieldIndex(i);
                const tensor D = dZones[j];
                const tensor F = fZones[j];

                const tensor Cd = mu[celli]*D + (rho[celli]*mag(U[celli]))*F;

                AU[celli] += Cd;

                porousGradp.primitiveFieldRef()[celli] = -(Cd &  U[celli]);
            }
        }

        porousGradp.correctBoundaryConditions();
    }
    else
    {
        forAll(cellZoneIDs_, zoneI)
        {
            const tensorField& dZones = D_[zoneI];
            const tensorField& fZones = F_[zoneI];

            const labelList& cells = mesh_.cellZones()[cellZoneIDs_[zoneI]];

            forAll(cells, i)
            {
                const label celli = cells[i];
                const label j = this->fieldIndex(i);
                const tensor D = dZones[j];
                const tensor F = fZones[j];

                AU[celli] += mu[celli]*D + (rho[celli]*mag(U[celli]))*F;
            }
        }
    }
}


// ************************************************************************* //
