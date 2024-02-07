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
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "farfieldRiemannFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Field<Type> farfieldRiemannFvPatchField<Type>::bcValue
(
            const fvPatchField<scalar>& pp,
            const fvPatchField<scalar>& prho,
            const fvPatchField<vector>& pU,
            const fvPatchField<scalar>& pT,
            const fluidThermo& thermo
) const
{
    const label patchi = this->patch().index();

    const tmp<volScalarField> tgamma = thermo.gamma();
    const scalarField& pgamma = tgamma().boundaryField()[patchi];

    const tmp<volScalarField> tpsi
    (
        this->db().objectRegistry::lookupObject<volScalarField>("thermo:psi")
    );
    const scalarField& ppsi = tpsi().boundaryField()[patchi];

    const scalarField R = 1.0/(ppsi*pT);

    const vectorField n = this->patch().nf();

    const scalarField pi = pp.patchInternalField();
    const scalarField rhoi = prho.patchInternalField();

    scalarField ci = sqrt(pgamma*pi/rhoi);
    scalarField si = pi*pow(1.0/rhoi, pgamma);
    const vectorField Ui = pU.patchInternalField();
    scalarField Mi = (Ui & n)/ci;

    scalarField cInf = sqrt(pgamma*R*TInf_);
    scalarField rhoInf = pgamma*pInf_/sqr(cInf);
    scalarField sInf = pInf_*pow(1.0/rhoInf, pgamma);
    vectorField UInf = MInf_*cInf*flowDir_;

    scalarField Rp = (Ui & n) + 2.0*ci/(pgamma - 1.0);
    scalarField Rm = (UInf & n) - 2.0*cInf/(pgamma - 1.0);

    forAll(Mi, facei) //for supersonic
    {
        if(Mi[facei] >= 1.0) //outflow
        {
            Rm[facei] = 
            (    
                (Ui[facei] & n[facei]) 
              - 2.0*ci[facei]/(pgamma[facei] - 1.0)
            );
        }
        else if(Mi[facei] <= -1.0) //inflow
        {
            Rp[facei] = 
            (
                (UInf[facei] & n[facei]) 
              + 2.0*cInf[facei]/(pgamma[facei] - 1.0)
            );
        }
    }

    scalarField cb = 0.25*(pgamma - 1.0)*(Rp - Rm);

    scalarField rhob(prho.size(), pTraits<scalar>::zero);
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
                    cb[facei]*cb[facei]/pgamma[facei]/si[facei], 
                    1.0/(pgamma[facei] - 1.0)
                );

            pb[facei] = rhob[facei]*sqr(cb[facei])/pgamma[facei];

            Tb[facei] = pb[facei]/rhob[facei]/R[facei];
        }
        else                //inflow
        {
            // using free-stream entropy
            rhob[facei] = 
                pow
                (
                    cb[facei]*cb[facei]/pgamma[facei]/sInf[facei], 
                    1.0/(pgamma[facei] - 1.0)
                );

            pb[facei] = rhob[facei]*sqr(cb[facei])/pgamma[facei];

            Tb[facei] = pb[facei]/rhob[facei]/R[facei];
        }
    }

    Field<Type> bcValue_(this->patch().size(), pTraits<Type>::zero);

    if(this->internalField().name() == "p")
    {
        bcValue_ = pb*pTraits<Type>::one;
    }
    else if(this->internalField().name() == "T")
    {
        bcValue_ = Tb*pTraits<Type>::one;
    }

    return bcValue_;
}

template<>
Field<vector> farfieldRiemannFvPatchField<vector>::bcValue
(
            const fvPatchField<scalar>& pp,
            const fvPatchField<scalar>& prho,
            const fvPatchField<vector>& pU,
            const fvPatchField<scalar>& pT,
            const fluidThermo& thermo
) const
{
    const label patchi = this->patch().index();

    const tmp<volScalarField> tgamma = thermo.gamma();
    const scalarField& pgamma = tgamma().boundaryField()[patchi];

    const tmp<volScalarField> tpsi
    (
        this->db().objectRegistry::lookupObject<volScalarField>("thermo:psi")
    );
    const scalarField& ppsi = tpsi().boundaryField()[patchi];

    const scalarField R = 1.0/(ppsi*pT);

    const vectorField n = this->patch().nf();

    const scalarField pi = pp.patchInternalField();
    const scalarField rhoi = prho.patchInternalField();

    scalarField ci = sqrt(pgamma*pi/rhoi);
    const vectorField Ui = pU.patchInternalField();
    scalarField Mi = (Ui & n)/ci;

    scalarField cInf = sqrt(pgamma*R*TInf_);
    vectorField UInf = MInf_*cInf*flowDir_;

    scalarField Rp = (Ui & n) + 2.0*ci/(pgamma - 1.0);
    scalarField Rm = (UInf & n) - 2.0*cInf/(pgamma - 1.0);

    forAll(Mi, facei) //for supersonic
    {
        if(Mi[facei] >= 1.0) //outflow
        {
            Rm[facei] = 
            (
                (Ui[facei] & n[facei]) 
              - 2.0*ci[facei]/(pgamma[facei] - 1.0)
            );
        }
        else if(Mi[facei] <= -1.0) //inflow
        {
            Rp[facei] = 
            (
                (UInf[facei] & n[facei])
              + 2.0*cInf[facei]/(pgamma[facei] - 1.0)
            );
        }
    }

    scalarField Un = 0.5*(Rp + Rm);

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

    Field<vector> bcValue_(this->patch().size(), pTraits<vector>::zero);

    if(this->internalField().name() == "U")
    {
        bcValue_ = Ub;
    }

    return bcValue_;
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
    flowDir_(Zero),
    curTimeIndex_(-1)
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
    MInf_(readScalar(dict.lookup("MInf"))),
    pInf_(readScalar(dict.lookup("pInf"))),
    TInf_(readScalar(dict.lookup("TInf"))),
    flowDir_(dict.lookup("flowDir")),
    curTimeIndex_(-1)
{
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
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    MInf_(ptf.MInf_),
    pInf_(ptf.pInf_),
    TInf_(ptf.TInf_),
    flowDir_(ptf.flowDir_),
    curTimeIndex_(-1)
{}


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
    flowDir_(ptf.flowDir_),
    curTimeIndex_(-1)
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
    flowDir_(ptf.flowDir_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void farfieldRiemannFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
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

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        Field<Type>& patchField = *this;

        const fvPatch& patch = this->patch();

        static word thermoName(fluidThermo::dictName);

        patchField = this->bcValue
        (
            patch.template lookupPatchField<volScalarField, scalar>("p"),
            patch.template lookupPatchField<volScalarField, scalar>("rho"),
            patch.template lookupPatchField<volVectorField, vector>("U"),
            patch.template lookupPatchField<volScalarField, scalar>("T"),
            this->db().objectRegistry::lookupObject<fluidThermo>(thermoName)
        );

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void farfieldRiemannFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("MInf") << MInf_ << token::END_STATEMENT << nl;
    os.writeKeyword("pInf") << pInf_ << token::END_STATEMENT << nl;
    os.writeKeyword("TInf") << TInf_ << token::END_STATEMENT << nl;
    os.writeKeyword("flowDir") << flowDir_ << token::END_STATEMENT << nl;
    this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
