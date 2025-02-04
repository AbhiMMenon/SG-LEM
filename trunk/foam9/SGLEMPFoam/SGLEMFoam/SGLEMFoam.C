/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2021 OpenFOAM Foundation
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
    Transient solver for turbulent flow of compressible reacting fluids with
    optional mesh motion and mesh topology changes.

    Uses the flexible PIMPLE (PISO-SIMPLE) solution for time-resolved and
    pseudo-transient simulations.

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
    #include "createNonPremixFields.H"
    #include "createSuperGrid.H"
    #include "internalFaces.H"
    #include "boundaryFaces.H"
    #include "initFiles.H"
    #include "initLEMlines.H"
    #include "initCmean.H"
    #include "initMatrices.H"
    #include "fillLEMlines.H"
    #include "integrateScalars.H" // initial mapping
    #include "createEnthalpy.H"

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
        #include "conInit.H"

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


        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
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

                if (pimple.firstPimpleIter() && !pimple.simpleRho())
                {
                    #include "rhoEqn.H"
                }

                if (pimple.models())
                {
                    fvModels.correct();
                }

                #include "UEqn.H"
                #include "ZEqn.H"
                #include "ZvarEqn.H"
                #include "CEqn.H"
                #include "CvarEqn.H"
                #include "integrateScalars.H"
        //      #include "computeHeatRelease.H"
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
                tStep = runTime.deltaTValue();
            }
        }

        rho = thermo.rho();

        Info << "\n\n-----------LEM advance + splicing-----------------\n\n";
        #include "wallSuperCells.H"
        #include "advanceLine.H"
//      #include "parallelChem.H"       // needs further testing, buggy
        #include "internalFlux.H"
        #include "boundaryFlux.H"
        #include "sorting.H"
        #include "spliceFrac.H"
        #include "inletSplice.H"
        #include "splice1.H"
        #include "parallelSplice.H"
        #include "splice2.H"
        Info << "\n--------------------------------------------------\n";

        runTime.write();
//      #include "writeLine.H"

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
