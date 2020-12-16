//##This script uses the rcnd command on the counting house machines to dump the most up to date info about all runs in the desired range. It is slow, but shouldn't need to be run often.

vector<TString> GetRunInfo(int run);

void MakeBigList(){
    
    //runrange to query. 2567-5000 prex, 5000-8572 crex
    int istart = 2566;
    int iend = 8573;
    vector<int> run_numbers;

    vector<vector<TString>> Table;
    for(int i = istart; i <= iend; i++){
        
        vector<TString> values = GetRunInfo(i);
        if(values[0] == "bad"){
            continue;
        }
        Table.push_back(values);
        run_numbers.push_back(i);
        
        //makes the status wiggle. partly for fun, partly as a visual cue.
        for(int j = 0; j<round(70*(sin((i-istart)*M_PI/20))+70); j++){
            cout << " ";
        }
        cout << "Finished run " << i << "..." << endl;
    }

    ofstream outputs;
    outputs.open("pcrex_run_data.list");
    
    //print table of values
    for(int i = 0; i < Table.size(); i++){
        int j = 0;
        for(j = 0; j < Table[i].size(); j++){
             if(j == Table[i].size()-1){
                 outputs << Table[i][j] << "," << run_numbers[i] << endl;
             }else{
                 outputs << Table[i][j] << ",";
             }
        }
    }
    outputs.close();
    
    return;
}

vector<TString> GetRunInfo(int run){
    
    //list of rcnd types, empty outputs default to space.
    vector<TString> types = {"arm_flag = ", "beam_current = ", "beam_energy = ", "bmw = ", "components = ", "component_stats = ", "event_count = ", "event_rate = ", "experiment = ", "feedback = ", "FFB = ", "flip_state = ", "good_charge = ", " helicity_frequency = ", "helicity_pattern = ", "horizontal_wien = ", "ihwp = ", "is_valid_run_end = ", "prompt_analysis = ", "respin_comment = ", "rhwp = ", "rtvs = ", "run_config = ", "run_end_time = ", "run_flag = ", "run_length = ", "run_prestart_time = ", "run_start_epoch = ", "run_start_time = ", "run_type = ", "session = ", "slug = ", "target_45encoder = ", "target_90encoder = ", "target_encoder = ", "target_type = ", "total_charge = ", "user_comment = ", "vertical_wien = ", "wac_comment = "};
    vector<TString> type_outputs = {" ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " ", " "};
    vector<TString> badrun = {"bad"};

    string output;
    output = gSystem->GetFromPipe(Form("rcnd %d",run));

    //turn long output into string stream to manipulate.
    stringstream outstream(output);
    while(1==1){
        string line;
        getline(outstream, line, '\n');
        
        if(line.find("is not found") != string::npos){
            return badrun;
        }
        
        //looks for the data type flags, splits the string at the first equals sign. Clears the output if it contains commas.
        for(int i = 0; i<types.size(); i++){
            if(line.find(types[i].Data()) != string::npos){
                int n = line.find("= ");
                string subline = line.substr(n+2, line.size());
                if(line.find(",") != string::npos){
                    subline = "BEGONE_COMMAS";
                }
                type_outputs[i] = subline;
            }
        }
        
        //break loop at end of stringstream. Add an end of line character becuse it output didnt have one.
        if(outstream.eof()){
            break;
        }
    }

    return type_outputs;
}
