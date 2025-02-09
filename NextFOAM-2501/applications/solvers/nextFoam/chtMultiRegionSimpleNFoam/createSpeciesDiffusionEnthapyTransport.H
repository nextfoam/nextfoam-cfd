// Calculate species diffusion entahlpy transport
// This functionality will need to be implemented later as part of a
// multiComponentMixture

if (Y.size())
{
    if (pvCoupling.solveEnergy() || pvCoupling.solveSpecies())
    {
        // This only applies to the dilute approximation for species mass 
        // diffusivities.This should be implemented later as a mixture model
        // library along with the unity Lewis number or multicomponent model.
        forAll(Y, i)
        {
            Dm.set
            (
                i, 
                new dimensionedScalar
                (
                    "Dm", 
                    dimArea/dimTime, 
                    1e-10, 
                    thermo.subDict(Y[i].name()).subDict("transport")
                )
            );   
        }

        Dt.reset(new volScalarField(turbulence->nut()/turbulence->Sct()));
    }

    if (pvCoupling.solveEnergy())
    {
        surfaceScalarField sumJ
        (
            surfaceScalarField::New
            (
                "sumJ",
                mesh,
                dimensionedScalar(dimMass/dimArea/dimTime, 0)
            )
        );

        surfaceScalarField sumJh
        (
            surfaceScalarField::New
            (
                "sumJh",
                mesh,
                dimensionedScalar(dimEnergy/dimArea/dimTime, 0)
            )
        );

        surfaceScalarField inertHsi
        (
            surfaceScalarField::New
            (
                "inertHsi",
                mesh,
                dimensionedScalar(dimEnergy/dimMass, 0)
            )
        );

        const volScalarField& T = thermo.T();

        forAll(Y, i)
        {
            const volScalarField& Yi = Y[i];

            tmp<volScalarField> tHsi
            (
                volScalarField::New
                (
                    IOobject::groupName("h", T.group()),
                    mesh,
                    dimEnergy/dimMass
                )
            );

            volScalarField& hsi = tHsi.ref();

            forAll(hsi, celli)
            {
                hsi[celli] = thermo.composition().Hs(i, p[celli], T[celli]);
            }

            volScalarField::Boundary& hsiBf = hsi.boundaryFieldRef();

            forAll(hsiBf, patchi)
            {
                const fvPatchScalarField& pT = T.boundaryField()[patchi];
                const fvPatchScalarField& pP = p.boundaryField()[patchi];

                forAll(hsiBf[patchi], patchFacei)
                {
                    hsiBf[patchi][patchFacei] = 
                        thermo.composition().Hs
                        (
                            i, 
                            pP[patchFacei], 
                            pT[patchFacei]
                        );
                }
            }

            if (i != inertIndex)
            {
                const surfaceScalarField ji
                (
                    fvc::interpolate(rho*(Dm[i] + Dt()))*fvc::snGrad(Yi)
                );

                sumJ += ji;

                sumJh += ji*fvc::interpolate(hsi);
            }
            else
            {
                inertHsi = fvc::interpolate(hsi);
            }
        }

        sumJh -= sumJ*inertHsi;

        speciesDiffusionEnthalpy.reset
        (
            new volScalarField(fvc::div(sumJh*mesh.magSf()))
        );
    }
}
