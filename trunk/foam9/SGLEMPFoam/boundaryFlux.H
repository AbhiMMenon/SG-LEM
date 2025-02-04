/*------------------------------------------------------------------------------
 * boundaryFluxes.H
 *
 * @author A.Menon 25th sep 2020
 *
 * MODIFIED 29th Sep
 *
 * Find boundaryFluxSC for all superGrid boundary faces; update the inFacesSC_loc
 * and outFacesSC_loc for all necessary cells at the boundary; make a
 * globalFlux list of both internal and boundary fluxes (loc numbering)
 *
 *-----------------------------------------------------------------------------*/

std::vector<double> boundaryFluxSC;
std::vector<double> boundaryCmeanSC;

for(int ix:patches){

    std::vector<double> Phi(phiVol.boundaryField()[ix].begin(),
                            phiVol.boundaryField()[ix].end());
    std::vector<double> Cbound(Cmean.boundaryField()[ix].begin(),
                            Cmean.boundaryField()[ix].end());
    
    for(size_t jx = 0; jx < boundaryOwnerSC[ix].size();jx++){

        double totFlux = 0.;
        double cLoc = 0.;
        int owner = boundaryOwnerSC[ix][jx];
        // add up individual CFD fluxes
        for(int facex:boundaryFaceListSC[ix][jx]){

            totFlux += Phi[facex];
            cLoc    += Cbound[facex];
        }

        cLoc /= boundaryFaceListSC[ix][jx].size();

        boundaryFluxSC.push_back(totFlux);
        boundaryCmeanSC.push_back(cLoc);
        int pos = boundaryFluxSC.size()-1;

        // assign superGrid boundary face number for influx or outflux
        if(totFlux < 0){
            //                            ↓  loc numbering  ↓
            inFacesSC_loc[owner].push_back(totFaces + pos); 
        }
        
        else {
            //                            ↓  loc numbering  ↓
           outFacesSC_loc[owner].push_back(totFaces + pos);
        }

    }

}

std::vector<double> globalFluxSC = fluxSC;
globalFluxSC.insert(globalFluxSC.end(), boundaryFluxSC.begin(), boundaryFluxSC.end()   );//← loc numbering

std::vector<double> globCfaceSC = CmeanSC_face;
globCfaceSC.insert(globCfaceSC.end(), boundaryCmeanSC.begin(),boundaryCmeanSC.end()   );//← loc numbering

Info << "SuperGrid boundary fluxes determined" << endl;
