/* Per wall patch, set wallPresent to 'true' for those LEM lines on the walls.
 * Find the super-grid averaged wall temperature and set lemLines[ix].wallT to
 * that temperature
 */

if(wallDiffusionFlag)
for(size_t ix = 0; ix < wallPatch.size(); ix++){
    int patchx = wallPatch[ix];

    // list the wall temperatures for this patch
    std::vector<double> Twall(T.boundaryField()[patchx].begin(),
                              T.boundaryField()[patchx].end());

    for(size_t jx = 0; jx < boundaryOwnerSC[patchx].size();jx++){
        int owner = boundaryOwnerSC[patchx][jx];

        // find average wall temperature for superCell
        double Tavg = 0.0;
        double count = 0.0;
        for(int facex:boundaryFaceListSC[patchx][jx]){
            Tavg += Twall[facex];
            count += 1.0;
        }
        Tavg /=count;
        lemLines[owner].wallPresent = true;
        lemLines[owner].wallT = Tavg;
        //test
        //lemLines[owner].setLEMT(Tavg);
    }

}
