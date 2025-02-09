/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2017 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "viscosityRatioInletOutletNuTildaFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::viscosityRatioInletOutletNuTildaFvPatchScalarField::
viscosityRatioInletOutletNuTildaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(p, iF),
    viscosityRatio_(nullptr)
{}


Foam::viscosityRatioInletOutletNuTildaFvPatchScalarField::
viscosityRatioInletOutletNuTildaFvPatchScalarField
(
    const viscosityRatioInletOutletNuTildaFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    inletOutletFvPatchScalarField(ptf, p, iF, mapper),
    viscosityRatio_(ptf.viscosityRatio_.clone(p.patch()))
{}


Foam::viscosityRatioInletOutletNuTildaFvPatchScalarField::
viscosityRatioInletOutletNuTildaFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    inletOutletFvPatchScalarField(p, iF),
    viscosityRatio_
    (
        PatchFunction1<scalar>::New(p.patch(), "viscosityRatio", dict)
    )
{
    fvPatchFieldBase::readDict(dict);
    this->phiName_ = dict.getOrDefault<word>("phi", "phi");   

    // The behavior when the value entry does not exist cannot be implemented
    // becaouse the decomposePar utility does not create a turbulence model 
    // obejct.
    this->readValueEntry(dict, IOobjectOption::MUST_READ);
}


Foam::viscosityRatioInletOutletNuTildaFvPatchScalarField::
viscosityRatioInletOutletNuTildaFvPatchScalarField
(
    const viscosityRatioInletOutletNuTildaFvPatchScalarField& ptf
)
:
    inletOutletFvPatchScalarField(ptf),
    viscosityRatio_(ptf.viscosityRatio_.clone(this->patch().patch()))
{}


Foam::viscosityRatioInletOutletNuTildaFvPatchScalarField::
viscosityRatioInletOutletNuTildaFvPatchScalarField
(
    const viscosityRatioInletOutletNuTildaFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    inletOutletFvPatchScalarField(ptf, iF),
    viscosityRatio_(ptf.viscosityRatio_.clone(this->patch().patch()))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::viscosityRatioInletOutletNuTildaFvPatchScalarField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    inletOutletFvPatchScalarField::autoMap(m);

    if (viscosityRatio_)
    {
        viscosityRatio_->autoMap(m);
    }
}


void Foam::viscosityRatioInletOutletNuTildaFvPatchScalarField::rmap
(
    const fvPatchScalarField& ptf,
    const labelList& addr
)
{
    inletOutletFvPatchScalarField::rmap(ptf, addr);

    const viscosityRatioInletOutletNuTildaFvPatchScalarField& tiptf =
        refCast<const viscosityRatioInletOutletNuTildaFvPatchScalarField>(ptf);

    if (viscosityRatio_)
    {
        viscosityRatio_->rmap(tiptf.viscosityRatio_(), addr);
    }
}


void Foam::viscosityRatioInletOutletNuTildaFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = db().time().timeOutputValue();

    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

    const fvsPatchScalarField& phip =
        patch().lookupPatchField<surfaceScalarField, scalar>(this->phiName_);

    this->refValue() = turbModel.nu(patch().index())*viscosityRatio_->value(t);
    this->valueFraction() = neg(phip);

    inletOutletFvPatchScalarField::updateCoeffs();
}


void Foam::viscosityRatioInletOutletNuTildaFvPatchScalarField::write(Ostream& os)
const
{
    fvPatchScalarField::write(os);
    
    if (viscosityRatio_)
    {
        viscosityRatio_->writeData(os);
    }

    os.writeEntryIfDifferent<word>("phi", "phi", this->phiName_);
    fvPatchField<scalar>::writeValueEntry(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        viscosityRatioInletOutletNuTildaFvPatchScalarField
    );
}

// ************************************************************************* //
