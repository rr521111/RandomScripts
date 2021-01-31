TGraphErrors* PlotFromMatrix(vector<vector<Double_t>> data, TString title);

void CompareTextOutputs(){

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

    ifstream infile;
    infile.open("./ComparisonOutputs/outputold.txt");
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

    ifstream infile2;
    infile2.open("./ComparisonOutputs/outputnew.txt");
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

    vector<vector<Double_t>> allruns;
    
    Double_t lastmini = 0;
    Double_t nextmini1 = 0;
    Double_t nextmini2 = 0;
    vector<Double_t> data1;
    vector<Double_t> data2;
    vector<Double_t> runmini;
    Int_t ministep = 5400;
    cout << maxcombo << endl;
    while(ministep < maxcombo+100){
        nextmini1 = 0;
        nextmini2 = 0;
        data1 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        data2 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
        runmini = {0, 0};
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

        ministep+=1;
        
        if(nextmini1 == maxcombo){
            break;
        }
    }
    
    Int_t totalcounts = 0;
    
    vector<vector<Double_t>> comparison;
    vector<vector<Double_t>> usedcomparison;
    vector<vector<Double_t>> asymold;
    vector<vector<Double_t>> asymnew;
    vector<vector<Double_t>> largecountdiffs;
    Int_t usedcountsold = 0;
    Int_t usedcountsnew = 0;
    vector<Int_t> sluglist = {100};
    cout << "Miniruns with absolute differences >= 13500:" << endl;
    for(Int_t i = 0; i < allruns.size(); i++){
        usedcountsold = 0;
        usedcountsnew = 0;
        comparison.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][3]-allruns[i][4], 0, 0});
        if(allruns[i][3] != 0 && allruns[i][5] != 0){
            asymold.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][5], 0, allruns[i][7]});
            usedcountsold = allruns[i][3];
        }
        if(allruns[i][4] != 0 && allruns[i][6] != 0){
            asymnew.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][6], 0, allruns[i][8]});
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

    //for(Int_t i = 0; i < sluglist.size(); i++){
    //    cout << sluglist[i] << endl;
    //}

    vector<vector<Double_t>> runcomparison;
    vector<vector<Double_t>> usedruncomparison;
    vector<vector<Double_t>> minicomparison;
    vector<vector<Double_t>> runasymold;
    vector<vector<Double_t>> runasymnew;
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
    cout << endl; 
    cout << "Runs with losses >= 10000:" << endl;
    for(Int_t i = 4000; i < 10000; i++){
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
            if(allruns[j][0] == i){
                numcountsold += allruns[j][3];
                numcountsnew += allruns[j][4];
                if(allruns[j][3] != 0){
                    if(allruns[j][5] != 0){
                        avgasymoldnum += allruns[j][5]/(allruns[j][7]*allruns[j][7]);
                        avgasymolddenom += 1/(allruns[j][7]*allruns[j][7]);
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

        if(numcountsold-numcountsnew >= 10000){
            cout << i << " " << numcountsold-numcountsnew << endl;
        }

        //cout << static_cast<Double_t>(i) << " " << avgasymoldnum/avgasymolddenom << " " << sqrt(1/avgasymolddenom) << endl;
        runcomparison.push_back({static_cast<Double_t>(i), static_cast<Double_t>(numcountsold-numcountsnew), 0, 0});
        usedruncomparison.push_back({static_cast<Double_t>(i), static_cast<Double_t>(usedcountsold-usedcountsnew), 0, 0});
        if(numminisold != 0 && usedcountsold != 0){
            runasymold.push_back({static_cast<Double_t>(i), avgasymoldnum/avgasymolddenom, 0, sqrt(1/avgasymolddenom)});
        }
        if(numminisnew != 0 && usedcountsnew != 0){
            runasymnew.push_back({static_cast<Double_t>(i), avgasymnewnum/avgasymnewdenom, 0, sqrt(1/avgasymnewdenom)});
        }
        minicomparison.push_back({static_cast<Double_t>(i), numminisold-numminisnew, 0, 0});
    }

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

        if(sluglist[i] <= 500){
            slugcomparison.push_back({static_cast<Double_t>(sluglist[i]), static_cast<Double_t>(numcountsold-numcountsnew), 0, 0});
            usedslugcomparison.push_back({static_cast<Double_t>(sluglist[i]), static_cast<Double_t>(usedcountsold-usedcountsnew), 0, 0});
            if(numminisold != 0 && usedcountsold != 0){
                slugasymold.push_back({static_cast<Double_t>(sluglist[i]), avgasymoldnum/avgasymolddenom, 0, sqrt(1/avgasymolddenom)});
            }
            if(numminisnew != 0 && usedcountsnew != 0){
                slugasymnew.push_back({static_cast<Double_t>(sluglist[i]), avgasymnewnum/avgasymnewdenom, 0, sqrt(1/avgasymnewdenom)});
            }
        }else{
            Double_t test = sluglist[i]%4000;
            slugcomparison.push_back({static_cast<Double_t>(sluglist[i])/15 + test, static_cast<Double_t>(numcountsold-numcountsnew), 0, 0});
            usedslugcomparison.push_back({static_cast<Double_t>(sluglist[i])/15 + test, static_cast<Double_t>(usedcountsold-usedcountsnew), 0, 0});
            if(numminisold != 0 && usedcountsold != 0){
                slugasymold.push_back({static_cast<Double_t>(sluglist[i])/15 + test, avgasymoldnum/avgasymolddenom, 0, sqrt(1/avgasymolddenom)});
            }
            if(numminisnew != 0 && usedcountsnew != 0){
                slugasymnew.push_back({static_cast<Double_t>(sluglist[i])/15 + test, avgasymnewnum/avgasymnewdenom, 0, sqrt(1/avgasymnewdenom)});
            }
        }
    }

    cout << "Total count difference: " << totalcounts1 << " " << totalcounts2 <<  " " << totalcounts1-totalcounts2 << " " << totalcounts << endl;
    
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
    TGraphErrors* asymoldplot = PlotFromMatrix(asymold, "Old Asymmetries Per Mini");
    asymoldplot->Draw("AP");
    asymoldplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c2->cd(2);
    TGraphErrors* asymnewplot = PlotFromMatrix(asymnew, "New Asymmetries Per Mini");
    asymnewplot->Draw("AP");
    asymnewplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c2->cd(3);
    TGraphErrors* runasymoldplot = PlotFromMatrix(runasymold, "Old Asymmetries Per Run");
    runasymoldplot->Draw("AP");
    runasymoldplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c2->cd(4);
    TGraphErrors* runasymnewplot = PlotFromMatrix(runasymnew, "New Asymmetries Per Run");
    runasymnewplot->Draw("AP");
    runasymnewplot->Fit("pol0");
    gStyle->SetOptFit(1); 
    c2->cd(5);
    TGraphErrors* slugasymoldplot = PlotFromMatrix(slugasymold, "Old Asymmetries Per Slug");
    slugasymoldplot->Draw("AP");
    slugasymoldplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c2->cd(6);
    TGraphErrors* slugasymnewplot = PlotFromMatrix(slugasymnew, "New Asymmetries Per Slug");
    slugasymnewplot->Draw("AP");
    slugasymnewplot->Fit("pol0");
    gStyle->SetOptFit(1);

    TCanvas *c3 = new TCanvas();
    TGraphErrors* minicountplot = PlotFromMatrix(minicomparison, "Minirun Count Differences Per Run");
    minicountplot->Draw("AP");
    minicountplot->Fit("pol0");
    gStyle->SetOptFit(1);
    
    /*cout << endl;
    cout << "Miniruns with differences >= 13500:" << endl;
    for(Int_t i = 0; i < largecountdiffs.size(); i++){
        cout << largecountdiffs[i][0] << " " << largecountdiffs[i][1] << endl;
    }*/

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
