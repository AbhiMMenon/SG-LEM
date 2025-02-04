/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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

Application
    reactingFoam

Description
    Solver for combustion with chemical reactions.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "fluidReactionThermo.H"
#include "combustionModel.H"
#include "dynamicMomentumTransportModel.H"
#include "fluidReactionThermophysicalTransportModel.H"
#include "multivariateScheme.H"
#include "pimpleControl.H"
#include "pressureReference.H"
#include "CorrectPhi.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "sampledSurfaces.H"
#include "cuttingPlane.H"

// LEM headers
#include "LEMLine.h"
#include "CellData.h"

// superGrid
#include "coarseMgMeshLevel.H"
#include "fineMgMeshLevel.H"

// splicing
#include "processorFvPatch.H"
#include "processorCyclicFvPatch.H"
#include "cyclicFvPatch.H"
#include "mpi.h"

#include <chrono>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{

    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "createDyMControls.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createRhoUfIfPresent.H"

    // <<-- SGLEM----------
    #include "createFields2.H"
    #include "createNumbers.H"
    #include "createSuperGrid.H"
    #include "internalFaces.H"
    #include "boundaryFaces.H"
    #include "initFiles.H"
    #include "initLEMlines.H"
    #include "initMatrices.H"
    #include "initCmean.H"
    #include "integrateScalars.H"
    thermo.correct();
    #include "createEnthalpy.H"
    #include "fillLEMlines.H"
    // ------------------->>

    turbulence->validate();

    if (!LTS)
    {
        #include "compressibleCourantNo.H"
        #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (pimple.run(runTime))
    {
        #include "readDyMControls.H"

        // Store divrhoU from the previous mesh so that it can be mapped
        // and used in correctPhi to ensure the corrected phi has the
        // same divergence
        autoPtr<volScalarField> divrhoU;
        if (correctPhi)
        {
            divrhoU = new volScalarField
            (
                "divrhoU",
                fvc::div(fvc::absolute(phi, rho, U))
            );
        }

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
            #include "compressibleCourantNo.H"
            #include "setDeltaT.H"
        }

        runTime++;
        
        double tStep ;
        tStep = runTime.deltaTValue();
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "computeNumbers.H"

        while (pimple.loop())
        {
            #include "conInit.H"
            if (!pimple.flow())
            {
                if (pimple.models())
                {
                    fvModels.correct();
                }

                if (pimple.thermophysics())
                {
                    #include "HEqn.H"
                }
            }
            else
            {
                if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
                {
                    // Store momentum to set rhoUf for introduced faces.
                    autoPtr<volVectorField> rhoU;
                    if (rhoUf.valid())
                    {
                        rhoU = new volVectorField("rhoU", rho*U);
                    }

                    fvModels.preUpdateMesh();

                    // Do any mesh changes
                    mesh.update();

                    if (mesh.changing())
                    {
                        MRF.update();

                        if (correctPhi)
                        {
                            #include "correctPhi.H"
                        }

                        if (checkMeshCourantNo)
                        {
                            #include "meshCourantNo.H"
                        }
                    }
                }

                if (pimple.firstPimpleIter() && !pimple.simpleRho())
                {
                    #include "rhoEqn.H"
                }

                if (pimple.models())
                {
                    fvModels.correct();
                }

                #include "UEqn.H"
                #include "CEqn.H"
                #include "CvarEqn.H"
                #include "integrateScalars.H"
                #include "HEqn.H"

                // --- Pressure corrector loop
                while (pimple.correct())
                {
                //  #include "../../compressible/rhoPimpleFoam/pEqn.H"
                    #include "pEqn.H"
                }

                if (pimple.turbCorr())
                {
                    turbulence->correct();
                    thermophysicalTransport->correct();
                }
            }
        }
        rho = thermo.rho();
        
        Info << "\n\n-----------LEM advance + splicing-----------------\n\n";
        #include "wallSuperCells.H"
        #include "advanceLine.H"
//      #include "parallelChem.H"
        #include "internalFlux.H"
        #include "boundaryFlux.H"
        #include "sorting.H"
        #include "spliceFrac.H"
        #include "inletSplice.H"
        #include "splice1.H"
        #include "parallelSplice.H"
        #include "splice2.H"
        Info << "\n--------------------------------------------------\n";
        
        bool wProps = runTime.write();
        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
   //   #include "writePlanes.H"
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
