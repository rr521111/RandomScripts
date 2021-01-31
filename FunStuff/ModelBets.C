vector<TString> types_line;
string ins;
vector<vector<Double_t>> games;

TGraphErrors* PlotFromMatrix(vector<vector<Double_t>> data, TString title);

void ModelBets() {
    TH1D* gamehist = new TH1D();
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

        cout << stod(types_line[1].Data()) << endl;
        if (types_line.size() <= 1) {
            break;
        }

        cout << stod(types_line[1].Data()) << endl;
        gamehist->Fill(stod(types_line[1].Data()));
        games.push_back({ stod(types_line[0].Data()), stod(types_line[1].Data()), 0, 0 });
    }

    TCanvas *c1 = new TCanvas();
    gamehist->Draw();

    TCanvas *c2 = new TCanvas();
    TGraphErrors *timegamehist = PlotFromMatrix(games, "First Million Crashes");
    timegamehist->Draw("AP");
    timegamehist->Fit("pol0");
    gStyle->SetOptFit(1);
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