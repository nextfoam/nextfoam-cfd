{
    const cellList& cells = mesh.cells();
    const label N = cells.size();
    const labelUList& own = mesh.owner();
    const labelUList& nei = mesh.neighbour();

   // Get the face area vector
    const surfaceVectorField& Sf = mesh.Sf();

    // calculate residual
    volSymmTensorField tau
    (
        "tau",
        -turbulence->devRhoReff()
        //-((2.0/3.0)*I)*rho*turbulence->k() //cwo
    );

    volScalarField rhoRes(fvc::div(Godunov->rhoFlux()));
    volVectorField rhoURes(fvc::div(Godunov->rhoUFlux()) - fvc::div(tau));
    volScalarField rhoERes
    (
        fvc::div(Godunov->rhoEFlux()) - fvc::div( tau & U )
      - fvc::laplacian(turbulence->alphaEff(), he)
      //- fvc::laplacian
      //  (
      //      (turbulence->mu() + 0.6*turbulence->mut()),
      //      turbulence->k()
      //  )
    );

    rhoRes.primitiveFieldRef() *= mesh.V();
    rhoURes.primitiveFieldRef() *= mesh.V();
    rhoERes.primitiveFieldRef() *= mesh.V();


    // get access to updated spectral radius
    const surfaceScalarField& spectralRadiiInv
    (
        spectralRadius.spectralRadiiInv()
    );

    const surfaceScalarField& spectralRadiiVis
    (
        spectralRadius.spectralRadiiVis()
    );

    const surfaceScalarField& spectralRadiiown
    (
        spectralRadius.spectralRadiiown()
    );

    const surfaceScalarField& spectralRadiinei
    (
        spectralRadius.spectralRadiinei()
    );

    // calcuate diagonal components
    scalarField diag(mesh.V()/deltaT);

    forAll(own, facei)
    {
        diag[own[facei]] += 
        (   
            0.5*spectralRadiiInv[facei] 
          + spectralRadiiVis[facei]/mesh.V()[own[facei]]
        );

        diag[nei[facei]] +=
        (
            0.5*spectralRadiiInv[facei] 
          + spectralRadiiVis[facei]/mesh.V()[nei[facei]]
        );
    }

    forAll(mesh.boundary(), patchi)
    {
        const labelUList& pFaceCells = mesh.boundary()[patchi].faceCells();

        const fvsPatchScalarField& pSpectralRadiiInv
        (
            spectralRadiiInv.boundaryField()[patchi]
        );

        const fvsPatchScalarField& pSpectralRadiiVis
        (
            spectralRadiiVis.boundaryField()[patchi]
        );

        forAll(mesh.boundary()[patchi], facei)
        {
            diag[pFaceCells[facei]] +=
            (
                0.5*pSpectralRadiiInv[facei] 
              + pSpectralRadiiVis[facei]/mesh.V()[pFaceCells[facei]]
            );
        }
    }

    scalarField& W1 = rho.primitiveFieldRef();
    vectorField& W2 = rhoU.primitiveFieldRef();
    scalarField& W3 = rhoE.primitiveFieldRef();
    //scalarField K = turbulence->k()().primitiveField(); //cwo

    const scalarField& W1_old = rho.oldTime().primitiveField();
    const vectorField& W2_old = rhoU.oldTime().primitiveField();
    const scalarField& W3_old = rhoE.oldTime().primitiveField();
    // turbulent kinetic energy will not updated in LU-SGS time integration
    //const scalarField& K_old = turbulence->k()().oldTime().primitiveField();
    //scalarField K_old = turbulence->k()().primitiveField(); //cwo

    const vectorField& cSf = Sf.primitiveField();

    // Forward sweep(Lower sweep)
    forAll(cells, celli)
    {
        const cell& cell = cells[celli];

        forAll(cell, facei)
        {
            label face = cell[facei];
            if(mesh.isInternalFace(face))
            {
                if(celli != own[face])
                {
                    label j = own[face];

                    scalar rhoFluxLower
                    (
                        ((W2_old[j] - W2[j]) & cSf[face])

                      - (W1[j] - W1_old[j])*spectralRadiiown[face]
                    );

                    vector rhoUFluxLower
                    (
                        (
                            (W2_old[j] & cSf[face])*W2_old[j]/W1_old[j]
                          + (
                                W3_old[j] 
                              - (W2_old[j] & W2_old[j])/2.0/W1_old[j]
                              //- W1_old[j]*K_old[j] //cwo
                            )*(gamma[j] - 1.0)*cSf[face]
                        )
                       
                      - (
                            (W2[j] & cSf[face])*W2[j]/W1[j]  
                          + (
                                W3[j] 
                              - (W2[j] & W2[j])/2.0/W1[j]
                              //- W1[j]*K[j] //cwo
                            )*(gamma[j] - 1.0)*cSf[face]
                        )

                      - (W2[j] - W2_old[j])*spectralRadiiown[face]
                    );

                    scalar rhoEFluxLower 
                    (
                        (
                            gamma[j]*W3_old[j]
                          - (
                                (W2_old[j] & W2_old[j])/2.0/W1_old[j]
                              //+ W1_old[j]*K_old[j] //cwo
                            )*(gamma[j] - 1.0)
                        )*(W2_old[j] & cSf[face])/W1_old[j]

                      - (
                            gamma[j]*W3[j]
                          - (
                                (W2[j] & W2[j])/2.0/W1[j]
                              //+ W1[j]*K[j] //cwo
                            )*(gamma[j] - 1.0)
                        )*(W2[j] & cSf[face])/W1[j]

                      - (W3[j] - W3_old[j])*spectralRadiiown[face]
                    );

                    W1[celli] -= (0.5*rhoFluxLower/diag[celli]); 
                    W2[celli] -= (0.5*rhoUFluxLower/diag[celli]); 
                    W3[celli] -= (0.5*rhoEFluxLower/diag[celli]); 
                }
            }
        }    

        W1[celli] -= (rhoRes.primitiveField()[celli]/diag[celli]);
        W2[celli] -= (rhoURes.primitiveField()[celli]/diag[celli]);
        W3[celli] -= (rhoERes.primitiveField()[celli]/diag[celli]);
    }

    // Backward sweep(Upper sweep)
    forAll(cells, i)
    {
        label celli = (N - 1) - i;
        const cell& c = cells[celli];

        forAll(c, facei)
        {
            label face = c[facei];
            if(mesh.isInternalFace(face))
            {
                if(celli != nei[face])
                {
                    label j = nei[face];

                    scalar rhoFluxUpper 
                    (
                        ((W2[j] - W2_old[j]) & cSf[face])
                      
                      - (W1[j] - W1_old[j])*spectralRadiinei[face]
                    );

                    vector rhoUFluxUpper
                    (
                        (        
                            (W2[j] & cSf[face])*W2[j]/W1[j]
                          + (
                                W3[j]  
                              - (W2[j] & W2[j])/2.0/W1[j]
                              //- W1[j]*K[j] //cwo
                            )*(gamma[j] - 1.0)*cSf[face]
                        )

                      - (
                            (W2_old[j] & cSf[face])*W2_old[j]/W1_old[j]
                          + (
                                W3_old[j]
                              - (W2_old[j] & W2_old[j])/2.0/W1_old[j]
                              //- W1_old[j]*K_old[j] //cwo
                            )*(gamma[j] - 1.0)*cSf[face]
                        )   
                                
                      - (W2[j] - W2_old[j])*spectralRadiinei[face]
                    );    

                    scalar rhoEFluxUpper
                    (
                        (
                            gamma[j]*W3[j]
                          - (
                                (W2[j] & W2[j])/2.0/W1[j]
                              //+ W1[j]*K[j] //cwo
                            )*(gamma[j] - 1.0)
                        )*(W2[j] & cSf[face])/W1[j]

                      - (
                            gamma[j]*W3_old[j]
                          - (
                                (W2_old[j] & W2_old[j])/2.0/W1_old[j]
                              //+ W1_old[j]*K_old[j] //cwo
                            )*(gamma[j] - 1.0)
                        )*(W2_old[j] & cSf[face])/W1_old[j]

                      - (W3[j] - W3_old[j])*spectralRadiinei[face]
                    );

                    W1[celli] -= 0.5*rhoFluxUpper/diag[celli];
                    W2[celli] -= 0.5*rhoUFluxUpper/diag[celli];
                    W3[celli] -= 0.5*rhoEFluxUpper/diag[celli];
                }
            }
        }
    }
}
