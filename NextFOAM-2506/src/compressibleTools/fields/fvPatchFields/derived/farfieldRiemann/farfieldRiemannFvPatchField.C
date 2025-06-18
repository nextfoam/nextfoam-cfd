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

\*---------------------------------------------------------------------------*/

#include "farfieldRiemannFvPatchField.H"

#include "thermodynamicConstants.H"
using namespace Foam::constant::thermodynamic;
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Field<Type> farfieldRiemannFvPatchField<Type>::bcValue
(
    const fvPatchField<scalar>& pp,
    const fvPatchField<vector>& pU,
    const fluidThermo& thermo
) const
{
    const label patchi = this->patch().index();

    const tmp<volScalarField> tgamma = thermo.gamma();
    const scalarField& gamma = tgamma().boundaryField()[patchi];

    const scalar molWeight
    (
        readScalar
        (
            thermo.subDict("mixture").subDict("specie").lookup("molWeight")
        )
    );

    const scalar R(RR/molWeight);

    const vectorField n(this->patch().nf());

    const scalarField pi(pp.patchInternalField());
    //const scalarField rhoi(prho.patchInternalField());
    const vectorField Ui(pU.patchInternalField());
    const scalarField Ti
    (
        thermo.T().boundaryField()[this->patch().index()].patchInternalField()
    );

    //scalarField ci(sqrt(gamma*pi/rhoi));
    scalarField ci(sqrt(gamma*R*Ti));
    //scalarField si(pi*pow(1.0/rhoi, gamma));
    scalarField si(pow(pi, 1.0 - gamma)*pow(R*Ti, gamma));
    scalarField Mi((Ui & n)/ci);

    scalarField cInf(sqrt(gamma*R*TInf_));
    scalarField rhoInf(gamma*pInf_/sqr(cInf));
    scalarField sInf(pInf_*pow(1.0/rhoInf, gamma));
    vectorField UInf(MInf_*cInf*flowDir_);

    scalarField Rp((Ui & n) + 2.0*ci/(gamma - 1.0));
    scalarField Rm((UInf & n) - 2.0*cInf/(gamma - 1.0));

    forAll(Mi, facei) //for supersonic
    {
        if(Mi[facei] >= 1.0) //outflow
        {
            Rm[facei] = 
            (    
                (Ui[facei] & n[facei]) 
              - 2.0*ci[facei]/(gamma[facei] - 1.0)
            );
        }
        else if(Mi[facei] <= -1.0) //inflow
        {
            Rp[facei] = 
            (
                (UInf[facei] & n[facei]) 
              + 2.0*cInf[facei]/(gamma[facei] - 1.0)
            );
        }
    }

    scalarField cb(0.25*(gamma - 1.0)*(Rp - Rm));

    scalarField rhob(pp.size(), pTraits<scalar>::zero);
    scalarField pb(pp.size(), pTraits<scalar>::zero);
    scalarField Tb(pp.size(), pTraits<scalar>::zero);

    forAll(Mi, facei)
    {
        if(Mi[facei] > 0.0) //outflow
        {
            // extrapolate entropy
            rhob[facei] = 
                pow
                (
                    cb[facei]*cb[facei]/gamma[facei]/si[facei], 
                    1.0/(gamma[facei] - 1.0)
                );

            pb[facei] = rhob[facei]*sqr(cb[facei])/gamma[facei];

            Tb[facei] = pb[facei]/rhob[facei]/R;
        }
        else                //inflow
        {
            // using free-stream entropy
            rhob[facei] = 
                pow
                (
                    cb[facei]*cb[facei]/gamma[facei]/sInf[facei], 
                    1.0/(gamma[facei] - 1.0)
                );

            pb[facei] = rhob[facei]*sqr(cb[facei])/gamma[facei];

            Tb[facei] = pb[facei]/rhob[facei]/R;
        }
    }

    Field<Type> bcValue(this->patch().size(), pTraits<Type>::zero);

    if
    (
        this->internalField().name() == "p"
     || this->internalField().name() == "p_rgh"
    )
    {
        bcValue = pb*pTraits<Type>::one;
    }
    else if(this->internalField().name() == "T")
    {
        bcValue = Tb*pTraits<Type>::one;
    }

    return bcValue;
}


template<>
Field<vector> farfieldRiemannFvPatchField<vector>::bcValue
(
    const fvPatchField<scalar>& pp,
    const fvPatchField<vector>& pU,
    const fluidThermo& thermo
) const
{
    const label patchi = this->patch().index();

    const tmp<volScalarField> tgamma = thermo.gamma();
    const scalarField& gamma = tgamma().boundaryField()[patchi];

    const scalar molWeight
    (
        readScalar
        (
            thermo.subDict("mixture").subDict("specie").lookup("molWeight")
        )
    );

    const scalar R(RR/molWeight);

    const vectorField n(this->patch().nf());

    const scalarField pi(pp.patchInternalField());
    //const scalarField rhoi(prho.patchInternalField());
    const vectorField Ui(pU.patchInternalField());
    const scalarField Ti
    (
        thermo.T().boundaryField()[this->patch().index()].patchInternalField()
    );

    //scalarField ci(sqrt(gamma*pi/rhoi));
    scalarField ci(sqrt(gamma*R*Ti));
    scalarField Mi((Ui & n)/ci);

    scalarField cInf(sqrt(gamma*R*TInf_));
    vectorField UInf(MInf_*cInf*flowDir_);

    scalarField Rp((Ui & n) + 2.0*ci/(gamma - 1.0));
    scalarField Rm((UInf & n) - 2.0*cInf/(gamma - 1.0));

    forAll(Mi, facei) //for supersonic
    {
        if(Mi[facei] >= 1.0) //outflow
        {
            Rm[facei] = 
            (
                (Ui[facei] & n[facei]) 
              - 2.0*ci[facei]/(gamma[facei] - 1.0)
            );
        }
        else if(Mi[facei] <= -1.0) //inflow
        {
            Rp[facei] = 
            (
                (UInf[facei] & n[facei])
              + 2.0*cInf[facei]/(gamma[facei] - 1.0)
            );
        }
    }

    scalarField Un(0.5*(Rp + Rm));

    vectorField Ub(pU.size(), pTraits<vector>::zero);

    forAll(Mi, facei)
    {
        if(Mi[facei] > 0.0) //outflow
        {
            Ub[facei] = 
                Ui[facei] + (Un[facei] - (Ui[facei] & n[facei]))*n[facei];
        }
        else //inflow
        {
            Ub[facei] = 
                UInf[facei] + (Un[facei] - (UInf[facei] & n[facei]))*n[facei];
        }
    }

    Field<vector> bcValue(this->patch().size(), pTraits<vector>::zero);

    if(this->internalField().name() == "U")
    {
        bcValue = Ub;
    }

    return bcValue;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
farfieldRiemannFvPatchField<Type>::farfieldRiemannFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    MInf_(0.0),
    pInf_(0.0),
    TInf_(0.0),
    flowDir_(Zero)
{}


template<class Type>
farfieldRiemannFvPatchField<Type>::farfieldRiemannFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const Field<Type>& fld
)
:
    fixedValueFvPatchField<Type>(p, iF, fld),
    MInf_(0.0),
    pInf_(0.0),
    TInf_(0.0),
    flowDir_(Zero)
{}


template<class Type>
farfieldRiemannFvPatchField<Type>::farfieldRiemannFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF, dict, false),
    MInf_(dict.get<scalar>("MInf")),
    pInf_(dict.get<scalar>("pInf")),
    TInf_(dict.get<scalar>("TInf")),
    flowDir_(dict.get<vector>("flowDir"))
{
    flowDir_ /= mag(flowDir_);

    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        this->evaluate();
    }
}


template<class Type>
farfieldRiemannFvPatchField<Type>::farfieldRiemannFvPatchField
(
    const farfieldRiemannFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(p, iF),    // Don't map
    MInf_(ptf.MInf_),
    pInf_(ptf.pInf_),
    TInf_(ptf.TInf_),
    flowDir_(ptf.flowDir_)
{
    if (mapper.direct() && !mapper.hasUnmapped())
    {
        // Use mapping instead of re-evaluation
        this->map(ptf, mapper);
    }
    else
    {
        // Evaluate since value not mapped
        this->evaluate();
    }
}


template<class Type>
farfieldRiemannFvPatchField<Type>::farfieldRiemannFvPatchField
(
    const farfieldRiemannFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    MInf_(ptf.MInf_),
    pInf_(ptf.pInf_),
    TInf_(ptf.TInf_),
    flowDir_(ptf.flowDir_)
{}


template<class Type>
farfieldRiemannFvPatchField<Type>::farfieldRiemannFvPatchField
(
    const farfieldRiemannFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    MInf_(ptf.MInf_),
    pInf_(ptf.pInf_),
    TInf_(ptf.TInf_),
    flowDir_(ptf.flowDir_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void farfieldRiemannFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& mapper
)
{
    fixedValueFvPatchField<Type>::autoMap(mapper);

    this->evaluate();
}


template<class Type>
void farfieldRiemannFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);
}


template<class Type>
void farfieldRiemannFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const fvPatch& patch = this->patch();

    static word thermoName(fluidThermo::dictName);

    fvPatchField<Type>::operator==
    (
        this->bcValue
        (
            patch.template lookupPatchField<volScalarField, scalar>("p"),
            patch.template lookupPatchField<volVectorField, vector>("U"),
            this->db().objectRegistry::template lookupObject<fluidThermo>
            (
                thermoName
            )
        )
    );

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void farfieldRiemannFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);

    os.writeEntry("MInf", MInf_);
    os.writeEntry("pInf", pInf_);
    os.writeEntry("TInf", TInf_);
    os.writeEntry("flowDir", flowDir_);

    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
