TGraphErrors* PlotFromMatrix(vector<vector<Double_t>> data, TString title);

void CompareTextOutputs(){

    vector<Int_t> runs1;
    vector<Int_t> miniruns1;
    vector<Int_t> counts1;
    vector<Double_t> asymmetries1;
    vector<Double_t> errors1;

    vector<Int_t> runs2;
    vector<Int_t> miniruns2;
    vector<Int_t> counts2;
    vector<Double_t> asymmetries2;
    vector<Double_t> errors2;

    Int_t totalcounts1 = 0;
    Int_t totalcounts2 = 0;

    Double_t maxcombo = 0;
    
    vector<TString> types_line = { " ", " ", " ", " ", " "};
    string ins;

    ifstream infile;
    infile.open("./outputold.txt");
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
        counts1.push_back(stoi(types_line[2].Data()));
        asymmetries1.push_back(stod(types_line[3].Data()));
        errors1.push_back(stod(types_line[4].Data()));
        //errors1.push_back(stod(types_line[4].Data())/sqrt(stod(types_line[2].Data())));
 
        totalcounts1+=stoi(types_line[2].Data());
        
        if(stod(types_line[0].Data())*100 + stod(types_line[1].Data()) > maxcombo){
            maxcombo = stod(types_line[0].Data())*100 + stod(types_line[1].Data());
        }
    }

    infile.close();

    ifstream infile2;
    infile2.open("./outputnew.txt");
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
        counts2.push_back(stoi(types_line[2].Data()));
        asymmetries2.push_back(stod(types_line[3].Data()));
        errors2.push_back(stod(types_line[4].Data()));
        //errors2.push_back(stod(types_line[4].Data())/sqrt(stod(types_line[2].Data())));
        
        totalcounts2+=stoi(types_line[2].Data());

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
        data1 = {0, 0, 0};
        data2 = {0, 0, 0};
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
                    data1 = {static_cast<Double_t>(counts1[i]), asymmetries1[i], errors1[i]};
                    runmini = {static_cast<Double_t>(runs1[i]), static_cast<Double_t>(miniruns1[i])};
                    break;
                }
            }
            for(Int_t i = 0; i < runs2.size(); i++){
                if(nextmini1==runs2[i]*100 + miniruns2[i]){
                    data2 = {static_cast<Double_t>(counts2[i]), asymmetries2[i], errors2[i]};
                    runmini = {static_cast<Double_t>(runs2[i]), static_cast<Double_t>(miniruns2[i])};
                    break;
                }
            }
            allruns.push_back({runmini[0], runmini[1], data1[0], data2[0], data1[1], data2[1], data1[2], data2[2]});
            
            lastmini=nextmini1;
        }else{
            for(Int_t i = 0; i < runs1.size(); i++){
                if(nextmini2==runs1[i]*100 + miniruns1[i]){
                    data1 = {static_cast<Double_t>(counts1[i]), asymmetries1[i], errors1[i]};
                    runmini = {static_cast<Double_t>(runs1[i]), static_cast<Double_t>(miniruns1[i])};
                    break;
                }
            }
            for(Int_t i = 0; i < runs2.size(); i++){
                if(nextmini2==runs2[i]*100 + miniruns2[i]){
                    data2 = {static_cast<Double_t>(counts2[i]), asymmetries2[i], errors2[i]};
                    runmini = {static_cast<Double_t>(runs2[i]), static_cast<Double_t>(miniruns2[i])};
                    break;
                }
            }
            allruns.push_back({runmini[0], runmini[1], data1[0], data2[0], data1[1], data2[1], data1[2], data2[2]});

            lastmini=nextmini2;
        }

        ministep+=1;
        
        if(nextmini1 == maxcombo){
            break;
        }
    }
    
    Int_t totalcounts = 0;
    
    vector<vector<Double_t>> comparison;
    vector<vector<Double_t>> asymold;
    vector<vector<Double_t>> asymnew;
    vector<vector<Double_t>> asymcomparison;
    for(Int_t i = 0; i < allruns.size(); i++){
        comparison.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][2]-allruns[i][3], 0, 0});
        asymold.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][4], 0, allruns[i][6]});
        asymnew.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][5], 0, allruns[i][7]});
        asymcomparison.push_back({allruns[i][0] + allruns[i][1]/20, allruns[i][4]-allruns[i][5], 0, sqrt(allruns[i][6]*allruns[i][6]+allruns[i][7]*allruns[i][7])});
        totalcounts+=allruns[i][2]-allruns[i][3];
    }

    cout << "Total count difference: " << totalcounts1 << " " << totalcounts2 <<  " " << totalcounts1-totalcounts2 << " " << totalcounts << endl;
    
    TCanvas *c1 = new TCanvas();
    TGraphErrors* countplot = PlotFromMatrix(comparison, "Count Differences");
    countplot->Draw("AP");
    countplot->Fit("pol0");
    gStyle->SetOptFit(1);
    
    TCanvas *c2 = new TCanvas();
    c2->Divide(2,1);
    c2->cd(1);
    TGraphErrors* asymoldplot = PlotFromMatrix(asymold, "Old Asymmetries");
    asymoldplot->Draw("AP");
    asymoldplot->Fit("pol0");
    gStyle->SetOptFit(1);
    c2->cd(2);
    TGraphErrors* asymnewplot = PlotFromMatrix(asymnew, "New Asymmetries");
    asymnewplot->Draw("AP");
    asymnewplot->Fit("pol0");
    gStyle->SetOptFit(1);

    TCanvas *c3 = new TCanvas();
    TGraphErrors* asymcomparisonplot = PlotFromMatrix(asymcomparison, "Asymmetry Differences");
    asymcomparisonplot->Draw("AP");
    
    
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
