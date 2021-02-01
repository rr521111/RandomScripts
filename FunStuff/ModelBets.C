vector<TString> types_line;
string ins;
vector<vector<Double_t>> games;

vector<vector<Double_t>> Model1(vector<vector<Double_t>> data, Double_t cutoff, Double_t betsize);

TGraphErrors* PlotFromMatrix(vector<vector<Double_t>> data, TString title);

void ModelBets() {
    Int_t i = 0;
    Double_t maxindex = -1;
    Double_t minindex = 99999999999999;
    Double_t maxmult = -1;
    Double_t minmult = 99999999999999;
    Double_t index = 0;
    Double_t mult = 0;
    Double_t random = 0;
    Double_t numcrash = 0;
    srand((unsigned) time(0));
    ifstream infile;
    infile.open("./crash_game_roobet.csv");
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

        if (types_line.size() <= 1) {
            break;
        }

        index = stod(types_line[0].Data());
        mult = stod(types_line[1].Data());

        if(i%100000==0){
            cout << "Read " << i << " events." << endl;
        }

        random = (rand() % 100) + 1;
        if(random <= 3){
            games.push_back({index, 0, 0, 0});
            numcrash++;
        }else{
            games.push_back({index, mult, 0, 0 });
        }

        if(mult>maxmult){
            maxmult=mult;
        }
    
        if(mult<minmult){
            minmult=mult;
        }

        if(index>maxindex){
            maxindex=index;
        }
        
        if(index<minindex){
            minindex=index;
        }
        i++;
    }
    cout << "Crash rate " << numcrash/i << endl;

    TH1D* gamehist = new TH1D("First Million Crashes", "First Million Crashes", 1000, minmult-1, maxmult+1);
    for(Int_t i = 0; i < games.size(); i++){
        gamehist->Fill(games[i][1]);
    }

    TCanvas *c1 = new TCanvas();
    gamehist->Draw();

    TH2D* gamehist2d = new TH2D("First Million Crashes Over Time", "First Million Crashes Over Time", 1000, minindex-1, maxindex+1, 1000, minmult-1, maxmult+1);
    for(Int_t i = 0; i < games.size(); i++){
        gamehist2d->Fill(games[i][0], games[i][1]);
    }

    TCanvas *c2 = new TCanvas();
    gamehist2d->Draw();
    gamehist2d->Fit("pol0");
    gStyle->SetOptFit(1);

    Double_t maxearn = -99999999999;
    Double_t minearn = 99999999999;
    vector<vector<Double_t>> Model1Earnings = Model1(games, 1.36, 1);
    for(Int_t i = 0; i < Model1Earnings.size(); i++){
        if(Model1Earnings[i][1]>maxearn){
            maxearn=Model1Earnings[i][1];
        }
        if(Model1Earnings[i][1]<minearn){
            minearn=Model1Earnings[i][1];
        }
    }

    TH2D* Model1hist2d = new TH2D("Model1 Earnings", "Model1 Earnings", 1000, minindex-1, maxindex+1, 1000, minearn-1, maxearn+1);
    for(Int_t i = 0; i < Model1Earnings.size(); i++){
        Model1hist2d->Fill(Model1Earnings[i][0], Model1Earnings[i][1]);
    }
     
    TCanvas *c3 = new TCanvas();
    Model1hist2d->Draw();
    Model1hist2d->Fit("pol1");
    gStyle->SetOptFit(1);
}

vector<vector<Double_t>> Model1(vector<vector<Double_t>> data, Double_t cutoff, Double_t betsize){
    vector<vector<Double_t>> earnings;

    for(Int_t i = 0; i < data.size(); i++){
        if(data[i][1]>cutoff){
            earnings.push_back({static_cast<Double_t>(i), data[i][1]*betsize, 0, 0});
        }else{
            earnings.push_back({static_cast<Double_t>(i), -1*betsize, 0, 0});
        }
    }

    return earnings;
};

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
