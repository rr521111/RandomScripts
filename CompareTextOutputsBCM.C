TGraphErrors* PlotFromMatrix(vector<vector<Double_t>> data, TString title);

void CompareTextOutputsBCM(){

    TString FileStub = "bcm";
    TString TargetString = "bcm";
    TString TreeVariable = "asym_bcm_an_us";

    vector<Int_t> runs1;
    vector<Int_t> miniruns1;
    vector<Int_t> slugs1;
    vector<Int_t> counts1;
    vector<Double_t> asymmetries1;
    vector<Double_t> errors1;
    vector<Double_t> bpm4ex1;
    vector<Double_t> bpm4exerrors1;
    vector<Double_t> bpm4ey1;
    vector<Double_t> bpm4eyerrors1;
    vector<Double_t> bpm4ax1;
    vector<Double_t> bpm4axerrors1;
    vector<Double_t> bpm4ay1;
    vector<Double_t> bpm4ayerrors1;

    vector<Int_t> runs2;
    vector<Int_t> miniruns2;
    vector<Int_t> slugs2;
    vector<Int_t> counts2;
    vector<Double_t> asymmetries2;
    vector<Double_t> errors2;
    vector<Double_t> bpm4ex2;
    vector<Double_t> bpm4exerrors2;
    vector<Double_t> bpm4ey2;
    vector<Double_t> bpm4eyerrors2;
    vector<Double_t> bpm4ax2;
    vector<Double_t> bpm4axerrors2;
    vector<Double_t> bpm4ay2;
    vector<Double_t> bpm4ayerrors2;

    Int_t totalcounts1 = 0;
    Int_t totalcounts2 = 0;

    Double_t maxcombo = 0;
    
    vector<TString> types_line = { " ", " ", " ", " ", " ", " "};
    string ins;

    //Reading in old data
    ifstream infile;
    infile.open(Form("./ComparisonOutputs/output%sold.txt", FileStub.Data()));
    while (!infile.eof()) {
        getline(infile, ins, '\n');
        stringstream line(ins);
        types_line.clear();
        //splits the line at commas
        while (line.good()) {
            string substr;
            getline(line, substr, ',');
            types_line.push_back(substr);
        }

        if(types_line.size() <= 1){
            break;
        }
        
        runs1.push_back(stoi(types_line[0].Data()));
        miniruns1.push_back(stoi(types_line[1].Data()));
        slugs1.push_back(stoi(types_line[2].Data()));
        counts1.push_back(stoi(types_line[3].Data()));
        asymmetries1.push_back(stod(types_line[4].Data()));
        errors1.push_back(stod(types_line[5].Data()));
        bpm4ex1.push_back(stod(types_line[6].Data()));
        bpm4exerrors1.push_back(stod(types_line[7].Data()));
        bpm4ey1.push_back(stod(types_line[8].Data()));
        bpm4eyerrors1.push_back(stod(types_line[9].Data()));
        bpm4ax1.push_back(stod(types_line[10].Data()));
        bpm4axerrors1.push_back(stod(types_line[11].Data()));
        bpm4ay1.push_back(stod(types_line[12].Data()));
        bpm4ayerrors1.push_back(stod(types_line[13].Data()));
        
 
        totalcounts1+=stoi(types_line[3].Data());
        
        if(stod(types_line[0].Data())*100 + stod(types_line[1].Data()) > maxcombo){
            maxcombo = stod(types_line[0].Data())*100 + stod(types_line[1].Data());
        }
    }

    infile.close();

    //Reading in new data
    ifstream infile2;
    infile2.open(Form("./ComparisonOutputs/output%snew.txt", FileStub.Data()));
    while (!infile2.eof()) {
        getline(infile2, ins, '\n');
        stringstream line(ins);
        types_line.clear();
        //splits the line at commas
        while (line.good()) {
            string substr;
            getline(line, substr, ',');
            types_line.push_back(substr);
        }

        if(types_line.size() <= 1){
            break;
        }
        
        runs2.push_back(stoi(types_line[0].Data()));
        miniruns2.push_back(stoi(types_line[1].Data()));
        slugs2.push_back(stoi(types_line[2].Data()));
        counts2.push_back(stoi(types_line[3].Data()));
        asymmetries2.push_back(stod(types_line[4].Data()));
        errors2.push_back(stod(types_line[5].Data()));
        bpm4ex2.push_back(stod(types_line[6].Data()));
        bpm4exerrors2.push_back(stod(types_line[7].Data()));
        bpm4ey2.push_back(stod(types_line[8].Data()));
        bpm4eyerrors2.push_back(stod(types_line[9].Data()));
        bpm4ax2.push_back(stod(types_line[10].Data()));
        bpm4axerrors2.push_back(stod(types_line[11].Data()));
        bpm4ay2.push_back(stod(types_line[12].Data()));
        bpm4ayerrors2.push_back(stod(types_line[13].Data()));
        
        totalcounts2+=stoi(types_line[3].Data());

        if(stod(types_line[0].Data())*100 + stod(types_line[1].Data()) > maxcombo){
            maxcombo = stod(types_line[0].Data())*100 + stod(types_line[1].Data());
        }
    }

    infile2.close();

    //Combining old and new data into one structure. Missing miniruns are zero'd.
    vector<vector<Double_t>> allruns;
    
    Double_t lastmini = 0;
    Double_t nextmini1 = 0;
    Double_t nextmini2 = 0;
    vector<Double_t> data1;
    vector<Double_t> data2;
    vector<Double_t> runmini;
    Int_t Step = 0;
    cout << maxcombo << endl;
    while(Step <= 100000){
        nextmini1 = 0;
        nextmini2 = 0;
        data1 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        data2 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        runmini = {0, 0};
        //Finds the next largest minirun between the two data sets. nextminis and combos shift the run number two decimals to have run/minirun in one number. 4000,11 becomes 400011
        for(Int_t i = 0; i < runs1.size(); i++){
            nextmini1 = runs1[i]*100 + miniruns1[i];
            if(nextmini1>lastmini){
                break;
            }
        }
        for(Int_t i = 0; i < runs2.size(); i++){
            nextmini2 = runs2[i]*100 + miniruns2[i];
            if(nextmini2>lastmini){
                break;
            }
        }

        //loops through both data sets to see if the next largest minirun matches data in either set. It should always find one match at least. If none are found for one dataset, it is left at 0.        
        if(nextmini1<nextmini2){
            for(Int_t i = 0; i < runs1.size(); i++){
                if(nextmini1==runs1[i]*100 + miniruns1[i]){
                    data1 = {static_cast<Double_t>(counts1[i]), asymmetries1[i], errors1[i], bpm4ex1[i], bpm4exerrors1[i], bpm4ey1[i], bpm4eyerrors1[i], bpm4ax1[i], bpm4axerrors1[i], bpm4ay1[i], bpm4ayerrors1[i]};
                    runmini = {static_cast<Double_t>(runs1[i]), static_cast<Double_t>(miniruns1[i]), static_cast<Double_t>(slugs1[i])};
                    break;
                }
            }
            for(Int_t i = 0; i < runs2.size(); i++){
                if(nextmini1==runs2[i]*100 + miniruns2[i]){
                    data2 = {static_cast<Double_t>(counts2[i]), asymmetries2[i], errors2[i], bpm4ex2[i], bpm4exerrors2[i], bpm4ey2[i], bpm4eyerrors2[i], bpm4ax2[i], bpm4axerrors2[i], bpm4ay2[i], bpm4ayerrors2[i]};
                    runmini = {static_cast<Double_t>(runs2[i]), static_cast<Double_t>(miniruns2[i]), static_cast<Double_t>(slugs2[i])};
                    break;
                }
            }
            
            allruns.push_back({runmini[0], runmini[1], runmini[2], data1[0], data2[0], data1[1], data2[1], data1[2], data2[2], data1[3], data2[3], data1[4], data2[4], data1[5], data2[5], data1[6], data2[6], data1[7], data2[7], data1[8], data2[8], data1[9], data2[9], data1[10], data2[10]});
            
            lastmini=nextmini1;

        }else{
            for(Int_t i = 0; i < runs1.size(); i++){
                if(nextmini2==runs1[i]*100 + miniruns1[i]){
                    data1 = {static_cast<Double_t>(counts1[i]), asymmetries1[i], errors1[i], bpm4ex1[i], bpm4exerrors1[i], bpm4ey1[i], bpm4eyerrors1[i], bpm4ax1[i], bpm4axerrors1[i], bpm4ay1[i], bpm4ayerrors1[i]};
                    runmini = {static_cast<Double_t>(runs1[i]), static_cast<Double_t>(miniruns1[i]), static_cast<Double_t>(slugs1[i])};
                    break;
                }
            }
            for(Int_t i = 0; i < runs2.size(); i++){
                if(nextmini2==runs2[i]*100 + miniruns2[i]){
                    data2 = {static_cast<Double_t>(counts2[i]), asymmetries2[i], errors2[i], bpm4ex2[i], bpm4exerrors2[i], bpm4ey2[i], bpm4eyerrors2[i], bpm4ax2[i], bpm4axerrors2[i], bpm4ay2[i], bpm4ayerrors2[i]};
                    runmini = {static_cast<Double_t>(runs2[i]), static_cast<Double_t>(miniruns2[i]), static_cast<Double_t>(slugs2[i])};
                    break;
                }
            }

            allruns.push_back({runmini[0], runmini[1], runmini[2], data1[0], data2[0], data1[1], data2[1], data1[2], data2[2], data1[3], data2[3], data1[4], data2[4], data1[5], data2[5], data1[6], data2[6], data1[7], data2[7], data1[8], data2[8], data1[9], data2[9], data1[10], data2[10]});

            lastmini=nextmini2;
        }

        //really complicated way to avoid missing or double counting the last minirun. big headache.
        if(Step > 0){
            if(allruns[Step-1][0]*100 + allruns[Step-1][1] == runmini[0]*100 + runmini[1]){
                cout << allruns[Step-1][0] << " " << allruns[Step-1][1] << " " << allruns[Step-1][2] << " " << allruns[Step-1][3] << " " << allruns[Step-1][4] << endl;
                cout << allruns[Step][0] << " " << allruns[Step][1] << " " << allruns[Step][2] << " " << allruns[Step][3] << " " << allruns[Step][4] << endl << endl;
                data1 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
                data2 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

                if(nextmini1>nextmini2){
                    nextmini2=nextmini1;
                }else{
                    nextmini1=nextmini2;
                }

                for(Int_t i = 0; i < runs1.size(); i++){
                    if(nextmini1==runs1[i]*100 + miniruns1[i]){
                        data1 = {static_cast<Double_t>(counts1[i]), asymmetries1[i], errors1[i], bpm4ex1[i], bpm4exerrors1[i], bpm4ey1[i], bpm4eyerrors1[i], bpm4ax1[i], bpm4axerrors1[i], bpm4ay1[i], bpm4ayerrors1[i]};
                        runmini = {static_cast<Double_t>(runs1[i]), static_cast<Double_t>(miniruns1[i]), static_cast<Double_t>(slugs1[i])};
                        break;
                    }
                }
                for(Int_t i = 0; i < runs2.size(); i++){
                    if(nextmini1==runs2[i]*100 + miniruns2[i]){
                        data2 = {static_cast<Double_t>(counts2[i]), asymmetries2[i], errors2[i], bpm4ex2[i], bpm4exerrors2[i], bpm4ey2[i], bpm4eyerrors2[i], bpm4ax2[i], bpm4axerrors2[i], bpm4ay2[i], bpm4ayerrors2[i]};
                        runmini = {static_cast<Double_t>(runs2[i]), static_cast<Double_t>(miniruns2[i]), static_cast<Double_t>(slugs2[i])};
                        break;
                    }
                }
            
                allruns[Step][0] = runmini[0];
                allruns[Step][1] = runmini[1];
                allruns[Step][2] = runmini[2];
                allruns[Step][3] = data1[0];
                allruns[Step][4] = data2[0];
                allruns[Step][5] = data1[1];
                allruns[Step][6] = data2[1];
                allruns[Step][7] = data1[2];
                allruns[Step][8] = data2[2];
                allruns[Step][9] = data1[3];
                allruns[Step][10] = data2[3];
                allruns[Step][11] = data1[4];
                allruns[Step][12] = data2[4];
                allruns[Step][13] = data1[5];
                allruns[Step][14] = data2[5];
                allruns[Step][15] = data1[6];
                allruns[Step][16] = data2[6];
                allruns[Step][17] = data1[7];
                allruns[Step][18] = data2[7];
                allruns[Step][19] = data1[8];
                allruns[Step][20] = data2[8];
                allruns[Step][21] = data1[9];
                allruns[Step][22] = data2[9];
                allruns[Step][23] = data1[10];
                allruns[Step][24] = data2[10];

                cout << allruns[Step-1][0] << " " << allruns[Step-1][1] << " " << allruns[Step-1][2] << " " << allruns[Step-1][3] << " " << allruns[Step-1][4] << endl;
                cout << allruns[Step][0] << " " << allruns[Step][1] << " " << allruns[Step][2] << " " << allruns[Step][3] << " " << allruns[Step][4] << endl << endl;

                break;
            }
        }

        Step += 1;
    }
    
    Int_t totalcounts = 0;
    
    //Making arrays for variables vs minirun.
    vector<vector<Double_t>> comparison;
    vector<vector<Double_t>> usedcomparison;
    vector<vector<Double_t>> asymold;
    vector<vector<Double_t>> asymnew;
    vector<vector<Double_t>> bpm4exnew;
    vector<vector<Double_t>> bpm4exold;
    vector<vector<Double_t>> bpm4eynew;
    vector<vector<Double_t>> bpm4eyold;
    vector<vector<Double_t>> bpm4axnew;
    vector<vector<Double_t>> bpm4axold;
    vector<vector<Double_t>> bpm4aynew;
    vector<vector<Double_t>> bpm4ayold;
    vector<vector<Double_t>> largecountdiffs;
    Int_t usedcountsold = 0;
    Int_t usedcountsnew = 0;
    vector<Int_t> sluglist = {100};
    cout << "Miniruns with absolute differences >= 13500:" << endl;
    for(Int_t i = 0; i < allruns.size(); i++){
        usedcountsold = 0;
        usedcountsnew = 0;
        comparison.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][3]-allruns[i][4], 0, 0});
        //Checks if count and asymmetry are nonzero.
        if(allruns[i][3] != 0 && allruns[i][5] != 0){
            asymold.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][5], 0, allruns[i][7]});
            bpm4exold.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][9]/1000000000, 0, allruns[i][11]/1000000000});
            bpm4eyold.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][13]/1000000000, 0, allruns[i][15]/1000000000});
            bpm4axold.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][17]/1000000000, 0, allruns[i][19]/1000000000});
            bpm4ayold.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][21]/1000000000, 0, allruns[i][23]/1000000000});
            usedcountsold = allruns[i][3];
        }
        if(allruns[i][4] != 0 && allruns[i][6] != 0){
            asymnew.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][6], 0, allruns[i][8]});
            bpm4exnew.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][10]/1000000000, 0, allruns[i][12]/1000000000});
            bpm4eynew.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][14]/1000000000, 0, allruns[i][16]/1000000000});
            bpm4axnew.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][18]/1000000000, 0, allruns[i][20]/1000000000});
            bpm4aynew.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][22]/1000000000, 0, allruns[i][24]/1000000000});
            usedcountsnew = allruns[i][4];
        }

        if(abs(usedcountsold-usedcountsnew) >= 13500){
            cout << allruns[i][0] << " " << allruns[i][1] << " " << usedcountsold-usedcountsnew << endl;
        }

        usedcomparison.push_back({allruns[i][0] + allruns[i][1]/20, static_cast<Double_t>(usedcountsold-usedcountsnew), 0, 0});
        totalcounts+=allruns[i][3]-allruns[i][4];
        for(Int_t j = 0; j < sluglist.size(); j++){
            if(allruns[i][2]  == sluglist[j]){
                break;
            }
            if(j == sluglist.size()-1){
                sluglist.push_back(allruns[i][2]);
            }
        }
    }

    //Arrays for variabls vs run.
    vector<vector<Double_t>> runcomparison;
    vector<vector<Double_t>> usedruncomparison;
    vector<vector<Double_t>> minicomparison;
    vector<vector<Double_t>> runasymold;
    vector<vector<Double_t>> runasymnew;
    vector<vector<Double_t>> runaqold;
    vector<vector<Double_t>> runaqnew;
    vector<vector<Double_t>> brokenminisold;
    vector<vector<Double_t>> brokenminisnew;
    Double_t numminisold = 0;
    Double_t numminisnew = 0;
    Int_t numcountsold = 0;
    Int_t numcountsnew = 0;
    Double_t avgasymoldnum = 0;
    Double_t avgasymolddenom = 0;
    Double_t avgasymnewnum = 0;
    Double_t avgasymnewdenom = 0;
    Double_t avgaqoldnum = 0;
    Double_t avgaqolddenom = 0;
    Double_t avgaqnewnum = 0;
    Double_t avgaqnewdenom = 0;
    cout << endl; 
    cout << "Runs with changes >= 10000:" << endl;
    for(Int_t i = 4000; i < 10000; i++){
        numcountsold = 0;
        numcountsnew = 0;
        usedcountsold = 0;
        usedcountsnew = 0;
        avgasymoldnum = 0;
        avgasymolddenom = 0;
        avgasymnewnum = 0;
        avgasymnewdenom = 0;
        avgaqoldnum = 0;
        avgaqolddenom = 0;
        avgaqnewnum = 0;
        avgaqnewdenom = 0;
        numminisold = 0;
        numminisnew = 0;
        for(Int_t j = 0; j < allruns.size(); j++){
            if(allruns[j][0] == i){
                numcountsold += allruns[j][3];
                numcountsnew += allruns[j][4];
                //Checks if counts and asym are nonzero, but counts the broken miniruns.
                if(allruns[j][3] != 0){
                    if(allruns[j][5] != 0){
                        avgasymoldnum += allruns[j][5]/(allruns[j][7]*allruns[j][7]);
                        avgasymolddenom += 1/(allruns[j][7]*allruns[j][7]);
                        avgaqoldnum += allruns[j][21]/(allruns[j][23]*allruns[j][23]);
                        avgaqolddenom += 1/(allruns[j][23]*allruns[j][23]);
                        usedcountsold += allruns[j][3];
                    }else{
                        brokenminisold.push_back(allruns[j]);
                    }
                    numminisold += 1;
                }
                if(allruns[j][4] != 0){
                    if(allruns[j][6] != 0){
                        avgasymnewnum += allruns[j][6]/(allruns[j][8]*allruns[j][8]);
                        avgasymnewdenom += 1/(allruns[j][8]*allruns[j][8]);
                        avgaqnewnum += allruns[j][22]/(allruns[j][24]*allruns[j][24]);
                        avgaqnewdenom += 1/(allruns[j][24]*allruns[j][24]);
                        usedcountsnew += allruns[j][4];
                    }else{
                        brokenminisnew.push_back(allruns[j]);
                    }
                    numminisnew += 1;
                }
            }
        }

        if(numminisold == 0 && numminisnew == 0){
            continue;
        }

        if(abs(numcountsold-numcountsnew) >= 10000){
            cout << i << " " << numcountsold-numcountsnew << endl;
        }

        runcomparison.push_back({static_cast<Double_t>(i), static_cast<Double_t>(numcountsold-numcountsnew), 0, 0});
        usedruncomparison.push_back({static_cast<Double_t>(i), static_cast<Double_t>(usedcountsold-usedcountsnew), 0, 0});
        if(numminisold != 0 && usedcountsold != 0){
            runasymold.push_back({static_cast<Double_t>(i), avgasymoldnum/avgasymolddenom, 0, sqrt(1/avgasymolddenom)});
            runaqold.push_back({static_cast<Double_t>(i), avgaqoldnum/avgaqolddenom, 0, sqrt(1/avgaqolddenom)});
        }
        if(numminisnew != 0 && usedcountsnew != 0){
            runasymnew.push_back({static_cast<Double_t>(i), avgasymnewnum/avgasymnewdenom, 0, sqrt(1/avgasymnewdenom)});
            runaqnew.push_back({static_cast<Double_t>(i), avgaqnewnum/avgaqnewdenom, 0, sqrt(1/avgaqnewdenom)});
        }
        minicomparison.push_back({static_cast<Double_t>(i), numminisold-numminisnew, 0, 0});
    }

    //Arrays for variables vs slug.
    vector<vector<Double_t>> slugcomparison;
    vector<vector<Double_t>> usedslugcomparison;
    vector<vector<Double_t>> slugasymold;
    vector<vector<Double_t>> slugasymnew;
    for(Int_t i = 0; i < sluglist.size(); i++){
        numcountsold = 0;
        numcountsnew = 0;
        usedcountsold = 0;
        usedcountsnew = 0;
        avgasymoldnum = 0;
        avgasymolddenom = 0;
        avgasymnewnum = 0;
        avgasymnewdenom = 0;
        numminisold = 0;
        numminisnew = 0;
        for(Int_t j = 0; j < allruns.size(); j++){
            if(allruns[j][2] == sluglist[i]){
                numcountsold += allruns[j][3];
                numcountsnew += allruns[j][4];
                if(allruns[j][3] != 0){
                    if(allruns[j][5] != 0){
                        avgasymoldnum += allruns[j][5]/(allruns[j][7]*allruns[j][7]);
                        avgasymolddenom += 1/(allruns[j][7]*allruns[j][7]);
                        usedcountsold += allruns[j][3];
                    }
                    numminisold += 1;
                }
                if(allruns[j][4] != 0){
                    if(allruns[j][6] != 0){
                        avgasymnewnum += allruns[j][6]/(allruns[j][8]*allruns[j][8]);
                        avgasymnewdenom += 1/(allruns[j][8]*allruns[j][8]);
                        usedcountsnew += allruns[j][4];
                    }
                    numminisnew += 1;
                }
            }
        }

        if(numminisold == 0 && numminisnew == 0){
            continue;
        }

        //if(sluglist[i] <= 500){
            slugcomparison.push_back({static_cast<Double_t>(sluglist[i]), static_cast<Double_t>(numcountsold-numcountsnew), 0, 0});
            usedslugcomparison.push_back({static_cast<Double_t>(sluglist[i]), static_cast<Double_t>(usedcountsold-usedcountsnew), 0, 0});
            if(numminisold != 0 && usedcountsold != 0){
                slugasymold.push_back({static_cast<Double_t>(sluglist[i]), avgasymoldnum/avgasymolddenom, 0, sqrt(1/avgasymolddenom)});
            }
            if(numminisnew != 0 && usedcountsnew != 0){
                slugasymnew.push_back({static_cast<Double_t>(sluglist[i]), avgasymnewnum/avgasymnewdenom, 0, sqrt(1/avgasymnewdenom)});
            }
        /*}else{
            Double_t test = sluglist[i]%4000;
            slugcomparison.push_back({static_cast<Double_t>(sluglist[i])/15 + test, static_cast<Double_t>(numcountsold-numcountsnew), 0, 0});
            usedslugcomparison.push_back({static_cast<Double_t>(sluglist[i])/15 + test, static_cast<Double_t>(usedcountsold-usedcountsnew), 0, 0});
            if(numminisold != 0 && usedcountsold != 0){
                slugasymold.push_back({static_cast<Double_t>(sluglist[i])/15 + test, avgasymoldnum/avgasymolddenom, 0, sqrt(1/avgasymolddenom)});
            }
            if(numminisnew != 0 && usedcountsnew != 0){
                slugasymnew.push_back({static_cast<Double_t>(sluglist[i])/15 + test, avgasymnewnum/avgasymnewdenom, 0, sqrt(1/avgasymnewdenom)});
            }
        }*/
    }

    cout << "Total count difference: " << totalcounts1 << " " << totalcounts2 <<  " " << totalcounts1-totalcounts2 << " " << totalcounts << endl;
    
    //Make a whole bunch of plots.
    TCanvas *c1 = new TCanvas();
    c1->Divide(1,3);
    c1->cd(1);
    TGraphErrors* countplot = PlotFromMatrix(comparison, "Count Differences Per Mini");
    countplot->Draw("AP");
    countplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c1->cd(2);
    TGraphErrors* runcountplot = PlotFromMatrix(runcomparison, "Count Differences Per Run");
    runcountplot->Draw("AP");
    runcountplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c1->cd(3);
    TGraphErrors* slugcountplot = PlotFromMatrix(slugcomparison, "Count Differences Per Slug");
    slugcountplot->Draw("AP");
    slugcountplot->Fit("pol0");
    gStyle->SetOptFit(1);

    TCanvas *c11 = new TCanvas();
    c11->Divide(1,3);
    c11->cd(1);
    TGraphErrors* countplot2 = PlotFromMatrix(usedcomparison, "Count Differences Per Mini No Broken Reg");
    countplot2->Draw("AP");
    countplot2->Fit("pol0");
    gStyle->SetOptFit(1);
    c11->cd(2);
    TGraphErrors* runcountplot2 = PlotFromMatrix(usedruncomparison, "Count Differences Per Run No Broken Reg");
    runcountplot2->Draw("AP");
    runcountplot2->Fit("pol0");
    gStyle->SetOptFit(1);
    c11->cd(3);
    TGraphErrors* slugcountplot2 = PlotFromMatrix(usedslugcomparison, "Count Differences Per Slug No Broken Reg");
    slugcountplot2->Draw("AP");
    slugcountplot2->Fit("pol0");
    gStyle->SetOptFit(1);
    
    TCanvas *c2 = new TCanvas();
    c2->Divide(2,3);
    c2->cd(1);
    TGraphErrors* asymoldplot = PlotFromMatrix(asymold, Form("Old %s %s Per Mini", TargetString.Data(), TreeVariable.Data()));
    asymoldplot->Draw("AP");
    asymoldplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c2->cd(2);
    TGraphErrors* asymnewplot = PlotFromMatrix(asymnew, Form("New %s %s Per Mini", TargetString.Data(), TreeVariable.Data()));
    asymnewplot->Draw("AP");
    asymnewplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c2->cd(3);
    TGraphErrors* runasymoldplot = PlotFromMatrix(runasymold, Form("Old %s %s Per Run", TargetString.Data(), TreeVariable.Data()));
    runasymoldplot->Draw("AP");
    runasymoldplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c2->cd(4);
    TGraphErrors* runasymnewplot = PlotFromMatrix(runasymnew, Form("New %s %s Per Run", TargetString.Data(), TreeVariable.Data()));
    runasymnewplot->Draw("AP");
    runasymnewplot->Fit("pol0");
    gStyle->SetOptFit(1); 
    c2->cd(5);
    TGraphErrors* slugasymoldplot = PlotFromMatrix(slugasymold, Form("Old %s %s Per Slug", TargetString.Data(), TreeVariable.Data()));
    slugasymoldplot->Draw("AP");
    slugasymoldplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c2->cd(6);
    TGraphErrors* slugasymnewplot = PlotFromMatrix(slugasymnew, Form("New %s %s Per Slug", TargetString.Data(), TreeVariable.Data()));
    slugasymnewplot->Draw("AP");
    slugasymnewplot->Fit("pol0");
    gStyle->SetOptFit(1);

    TCanvas *c3 = new TCanvas();
    TGraphErrors* minicountplot = PlotFromMatrix(minicomparison, "Minirun Count Differences Per Run");
    minicountplot->Draw("AP");
    minicountplot->Fit("pol0");
    gStyle->SetOptFit(1);

    TCanvas *c4 = new TCanvas();
    c4->Divide(2,4);
    c4->cd(1);
    TGraphErrors* bpm4exoldplot = PlotFromMatrix(bpm4exold, "Old bpm4eX yield Per Mini");
    bpm4exoldplot->Draw("AP");
    bpm4exoldplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c4->cd(2);
    TGraphErrors* bpm4exnewplot = PlotFromMatrix(bpm4exnew, "New bpm4eX yield Per Mini");
    bpm4exnewplot->Draw("AP");
    bpm4exnewplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c4->cd(3);
    TGraphErrors* bpm4eyoldplot = PlotFromMatrix(bpm4eyold, "Old bpm4eY yield Per Mini");
    bpm4eyoldplot->Draw("AP");
    bpm4eyoldplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c4->cd(4);
    TGraphErrors* bpm4eynewplot = PlotFromMatrix(bpm4eynew, "New bpm4eY yield Per Mini");
    bpm4eynewplot->Draw("AP");
    bpm4eynewplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c4->cd(5);
    TGraphErrors* bpm4axoldplot = PlotFromMatrix(bpm4axold, "Old bpm4aX yield Per Mini");
    bpm4axoldplot->Draw("AP");
    bpm4axoldplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c4->cd(6);
    TGraphErrors* bpm4axnewplot = PlotFromMatrix(bpm4axnew, "New bpm4aX yield Per Mini");
    bpm4axnewplot->Draw("AP");
    bpm4axnewplot->Fit("pol0");
    gStyle->SetOptFit(1);

    TCanvas *c5 = new TCanvas();
    c5->Divide(2,2);
    c5->cd(1);
    TGraphErrors* bpm4ayoldplot = PlotFromMatrix(bpm4ayold, "Old Aq Per Mini");
    bpm4ayoldplot->Draw("AP");
    bpm4ayoldplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c5->cd(2);
    TGraphErrors* bpm4aynewplot = PlotFromMatrix(bpm4aynew, "New Aq Per Mini");
    bpm4aynewplot->Draw("AP");
    bpm4aynewplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c5->cd(3);
    TGraphErrors* runaqoldplot = PlotFromMatrix(runaqold, "Old Aq Per Run");
    runaqoldplot->Draw("AP");
    runaqoldplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c5->cd(4);
    TGraphErrors* runaqnewplot = PlotFromMatrix(runaqnew, "New Aq Per Run");
    runaqnewplot->Draw("AP");
    runaqnewplot->Fit("pol0");
    gStyle->SetOptFit(1);

    cout << endl;
    cout << "Old miniruns with failed reg:" << endl;
    for(Int_t i = 0; i < brokenminisold.size(); i++){
        cout << brokenminisold[i][0] << " " << brokenminisold[i][1] << endl;
    }

    cout << endl;
    cout << "New miniruns with failed reg:" << endl;
    for(Int_t i = 0; i < brokenminisnew.size(); i++){
        cout << brokenminisnew[i][0] << " " << brokenminisnew[i][1] << endl;
    }
    
}

TGraphErrors* PlotFromMatrix(vector<vector<Double_t>> data, TString title) {

    Double_t x[data.size()];
    Double_t y[data.size()];
    Double_t ex[data.size()];
    Double_t ey[data.size()];
    for (Int_t i = 0; i < data.size(); i++) {
        x[i] = data[i][0];
        y[i] = data[i][1];
        ex[i] = data[i][2];
        ey[i] = data[i][3];
    }

    TGraphErrors* output = new TGraphErrors(data.size(), x, y, ex, ey);
    output->SetTitle(title);
    output->SetMarkerStyle(7);

    return output;
}
