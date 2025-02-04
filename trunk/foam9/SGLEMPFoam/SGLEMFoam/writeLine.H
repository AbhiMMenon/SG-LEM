/* Hard coded for MRB_2D case for diagnosing extinction
 *
 */

if(Pstream::myProcNo() ==8){
    int outL = 40;
    if(runTime.timeIndex()%10==0){
  //if(true){
        std::string FILENAME = "./LineProps/" + std::to_string(fileNo) + ".dat";
        file.open(FILENAME);
        lemLines[outL].writeProps();
        file.close();
        fileNo++;
    }
    
    if(runTime.timeIndex()%100==0){
        std::string FILENAME = "./LinePropsZC/" + std::to_string(fileNo) + ".dat";
        file.open(FILENAME);
        lemLines[outL].writeZCProps();
        file.close();
     // fileNo++;
    }

}
