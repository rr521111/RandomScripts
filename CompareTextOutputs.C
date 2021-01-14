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
        
        runs1.push_back(stod(types_line[0].Data()));
        counts1.push_back(stoi(types_line[1].Data()));
        asymmetries1.push_back(stod(types_line[2].Data()));
        errors1.push_back(stod(types_line[3].Data()));
        
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
        
        runs2.push_back(stod(types_line[0].Data()));
        counts2.push_back(stoi(types_line[1].Data()));
        asymmetries2.push_back(stod(types_line[2].Data()));
        errors2.push_back(stod(types_line[3].Data()));
        
    }

    infile2.close();

    Int_t longer;
    if(runs1.size() >= runs1.size()){
        longer = runs1.size();
    }else{
        longer = runs2.size();
    }

    Int_t totalcounts = 0;
    Int_t missing1 = 0;
    Int_t missing2 = 0;
    Int_t mulspermini = 9000;
    vector<vector<Double_t>> comparison;
    for(Int_t i = 0; i<longer; i++){
        if(i+missing1>runs1.size() || i+missing2>runs2.size()){
            break;
        }

        if(runs1[i+missing1]!=runs2[i+missing2]){
            cout << runs1[i+missing1] << " " << runs2[i+missing2] << endl;
        }

        if(runs1[i+missing1]==runs2[i+missing2]){
            comparison.push_back({runs1[i+missing1],static_cast<Double_t>(counts1[i+missing1]-counts2[i+missing2]),0,0});
            totalcounts += counts1[i+missing1]-counts2[i+missing2];
            cout << runs1[i+missing1] << " " << runs2[i+missing2] << "   " << static_cast<Double_t>(counts1[i+missing1]-counts2[i+missing2]) << endl;
        }else if(runs1[i+missing1] > runs2[i+missing2]){
            missing2++;
            i--;
        }else if(runs1[i+missing1] < runs2[i+missing2]){
            missing1++;
            i--;
        }
    }

    cout << "Total count difference: " << totalcounts << endl;
    
    TGraphErrors* countplot = PlotFromMatrix(comparison, "Count Differences");

    countplot->Draw("AP");
    
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
