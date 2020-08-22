vector<vector<Double_t>> FindNearestNeighbors(vector<Double_t> vertex, vector<vector<Double_t>> current_solid);

Bool_t IsFace(vector<vector<Double_t>> triplet, vector<vector<Double_t>> current_solid);

vector<vector<Double_t>> NormalizePoints(vector<vector<Double_t>> current_solid);

vector<vector<Double_t>> Resample(vector<vector<Double_t>> current_solid, vector<vector<Int_t>> current_edges);

vector<vector<Int_t>> FindMissingEdges(Int_t oldverts, vector<vector<Double_t>> current_solid, vector<vector<Int_t>> current_edges);

Int_t GetVertexNumber(vector<Double_t> vertex, vector<vector<Double_t>> current_solid);

vector<Double_t> ConvertToCart(vector<Double_t> vector);

vector<Double_t> ConvertToSpherical(vector<Double_t> vector);

Double_t r = 1000;
Int_t Iterations = 2; 

void SphereMapping(){
    vector<vector<Double_t>> Icosahedron = {
        {0, 1, (1+sqrt(5))/sqrt(2)},
        {0, -1, (1+sqrt(5))/sqrt(2)},
        {(1+sqrt(5))/sqrt(2), 0, 1},
        {-(1+sqrt(5))/sqrt(2), 0, 1},
        {1, (1+sqrt(5))/sqrt(2), 0},
        {-1, (1+sqrt(5))/sqrt(2), 0},
        {1, -(1+sqrt(5))/sqrt(2), 0},
        {-1, -(1+sqrt(5))/sqrt(2), 0},
        {(1+sqrt(5))/sqrt(2), 0, -1},
        {-(1+sqrt(5))/sqrt(2), 0, -1},
        {0, 1, -(1+sqrt(5))/sqrt(2)},
        {0, -1, -(1+sqrt(5))/sqrt(2)}
    };
    
    vector<vector<Int_t>> Icosahedron_Edges = {
        {0, 1},
        {1, 2},
        {2, 0},
        {0, 3},
        {3, 1},
        {1, 6},
        {6, 7},
        {7, 1},
        {0, 4},
        {4, 5},
        {5, 0},
        {4, 2},
        {2, 6},
        {6, 11},
        {11, 7},
        {7, 3},
        {3, 9},
        {9, 7},
        {9, 5},
        {5, 3},
        {5, 10},
        {10, 4},
        {4, 8},
        {8, 2},
        {8, 6},
        {11, 10},
        {10, 9},
        {9, 11},
        {11, 8},
        {8, 10}
    };
    
    Icosahedron = NormalizePoints(Icosahedron);
    cout << "Verticies: " << Icosahedron.size() << ",  Edges: " << Icosahedron_Edges.size() << endl;

    Int_t verticies = Icosahedron.size();
    for(Int_t i = 0; i < Iterations; i++){
        Icosahedron = Resample(Icosahedron, Icosahedron_Edges);
        Icosahedron = NormalizePoints(Icosahedron);
        Icosahedron_Edges = FindMissingEdges(verticies, Icosahedron, Icosahedron_Edges);
        cout << "Verticies: " << Icosahedron.size() << ",  Edges: " << Icosahedron_Edges.size() << endl;
        verticies = Icosahedron.size();
    }

    vector<TPolyLine3D*> test;
    for(Int_t i = 0; i < Icosahedron_Edges.size(); i++){
        test.push_back(new TPolyLine3D());
        test[i]->SetPoint(0, Icosahedron[Icosahedron_Edges[i][0]][0], Icosahedron[Icosahedron_Edges[i][0]][1], Icosahedron[Icosahedron_Edges[i][0]][2]);
        test[i]->SetPoint(1, Icosahedron[Icosahedron_Edges[i][1]][0], Icosahedron[Icosahedron_Edges[i][1]][1], Icosahedron[Icosahedron_Edges[i][1]][2]);
    }

    TH3D* visual = new TH3D("", "", 100, -1000, 1000, 100, -1000, 1000, 100, -1000, 1000);

    for(Int_t i = 0; i < Icosahedron.size(); i++){
        visual->Fill(Icosahedron[i][0], Icosahedron[i][1], Icosahedron[i][2]);
    }

    gStyle->SetCanvasPreferGL(1);
    visual->SetMarkerStyle(7);

    test[0]->Draw();
    for(Int_t i = 1; i < test.size(); i++){
        test[i]->Draw("same");
    }
    //visual->Draw();
}

vector<vector<Double_t>> FindNearestNeighbors(vector<Double_t> vertex, vector<vector<Double_t>> current_solid){
    vector<vector<Double_t>> Neighbors;
   
    Double_t distance;
    Double_t min_distance = 10000;
    
    for(Int_t i = 0; i < current_solid.size(); i++){
        distance = sqrt(pow((vertex[0]-current_solid[i][0]), 2) + pow((vertex[1]-current_solid[i][1]), 2) + pow((vertex[2]-current_solid[i][2]), 2));
        if(distance < 10){
            continue;
        }else if(distance < min_distance){
            min_distance = distance;
        }
    }
    
    for(Int_t radius = 0; radius < r; radius++){
        for(Int_t i = 0; i < current_solid.size(); i++){
            distance = sqrt(pow((vertex[0]-current_solid[i][0]), 2) + pow((vertex[1]-current_solid[i][1]), 2) + pow((vertex[2]-current_solid[i][2]), 2));
            if(distance <= min_distance+radius && distance > 1){
                Neighbors.push_back(current_solid[i]);
            }
        }

        if(Neighbors.size()==6){
           break;
        }

        Neighbors.clear();
    }
    
    return Neighbors;
}

Bool_t IsFace(vector<vector<Double_t>> triplet, vector<vector<Double_t>> current_solid){
    
    Double_t distance;
    Double_t min_distance;
    
    for(Int_t i = 0; i < current_solid.size(); i++){
        distance = sqrt(pow((triplet[0][0]-current_solid[i][0]), 2) + pow((triplet[0][1]-current_solid[i][1]), 2) + pow((triplet[0][2]-current_solid[i][2]), 2));
        if(distance < min_distance){
            min_distance = distance;
        }
    }

    Double_t AB = sqrt(pow((triplet[0][0]-triplet[1][0]), 2) + pow((triplet[0][1]-triplet[1][1]), 2) + pow((triplet[0][2]-triplet[1][2]), 2));
    Double_t BC = sqrt(pow((triplet[2][0]-triplet[1][0]), 2) + pow((triplet[2][1]-triplet[1][1]), 2) + pow((triplet[2][2]-triplet[1][2]), 2));
    Double_t CA = sqrt(pow((triplet[0][0]-triplet[2][0]), 2) + pow((triplet[0][1]-triplet[2][1]), 2) + pow((triplet[0][2]-triplet[2][2]), 2));
    
    if(AB == BC && BC == CA){
        return true;
    }else{
        return false;
    }
}

vector<vector<Double_t>> NormalizePoints(vector<vector<Double_t>> current_solid){
    vector<vector<Double_t>> output;
    vector<Double_t> cart;
    vector<Double_t> sphere;

    for(Int_t i = 0; i < current_solid.size(); i++){
        sphere = ConvertToSpherical(current_solid[i]);
        sphere[0] = r;
        cart = ConvertToCart(sphere);
        output.push_back(cart);
    }

    return output;
}

vector<vector<Double_t>> Resample(vector<vector<Double_t>> current_solid, vector<vector<Int_t>> current_edges){
    Double_t x;
    Double_t y;
    Double_t z;

    for(Int_t i = 0; i < current_edges.size(); i++){ 
        x = (current_solid[current_edges[i][0]][0] + current_solid[current_edges[i][1]][0])/2;
        y = (current_solid[current_edges[i][0]][1] + current_solid[current_edges[i][1]][1])/2;
        z = (current_solid[current_edges[i][0]][2] + current_solid[current_edges[i][1]][2])/2;

        current_solid.push_back({x,y,z});
    }

    return current_solid;
}

vector<vector<Int_t>> FindMissingEdges(Int_t oldverts, vector<vector<Double_t>> current_solid, vector<vector<Int_t>> current_edges){
    vector<vector<Int_t>> new_edges;
    
    for(Int_t i = 0; i < current_edges.size(); i++){ 
        new_edges.push_back({current_edges[i][0], i + oldverts});
        new_edges.push_back({current_edges[i][1], i + oldverts});
    }

    vector<vector<Int_t>> new_edges2;
    vector<vector<Double_t>> NearestNeighbors;
    for(Int_t i = 0; i < current_solid.size(); i++){
        NearestNeighbors = FindNearestNeighbors(current_solid[i], current_solid);
        cout << NearestNeighbors.size() << endl;
        
        for(Int_t k = 0; k < NearestNeighbors.size(); k++){
            Bool_t exists = false;
            vector<Int_t> test_edge1 = {GetVertexNumber(NearestNeighbors[k], current_solid), i};
            vector<Int_t> test_edge2 = {i, GetVertexNumber(NearestNeighbors[k], current_solid)};

            for(Int_t j = 0; j < new_edges.size(); j++){
                if(test_edge1 == new_edges[j] || test_edge2 == new_edges[j]){
                    exists = true;
                }
            }
      
            if(!exists){
                new_edges2.push_back(test_edge1);
            }
        }

        if(current_solid.size()>42){
            if(current_solid.size()>168){
                if(current_solid.size()>696){
                    Int_t one_percent = round(current_solid.size()/100);
                    if(0 == i%one_percent){
                        cout << "About " << round(i/one_percent) << "% done!" << endl;
                    }
                    continue;
                }

                Int_t five_percent = round(current_solid.size()/20);
                if(0 == i%five_percent){
                    cout << "About " << 5*round(i/five_percent) << "% done!" << endl;
                }
                continue;
            }
            
            Int_t ten_percent = round(current_solid.size()/10);
            if(0 == i%ten_percent){
                cout << "About " << 10*round(i/ten_percent) << "% done!" << endl;
            }
        }
    }
 
    for(Int_t i = 0; i < new_edges2.size(); i++){
        Bool_t exists = false;
        vector<Int_t> temp = {new_edges2[i][1], new_edges2[i][0]};

        for(Int_t j = 0; j < new_edges.size(); j++){
            if(new_edges2[i] == new_edges[j] || temp == new_edges[j]){
                exists = true;
            }
        }
    
        if(!exists){
            new_edges.push_back(new_edges2[i]);
        }
    }
    
    return new_edges;
}

Int_t GetVertexNumber(vector<Double_t> vertex, vector<vector<Double_t>> current_solid){
    Int_t number;
    for(Int_t i = 0; i < current_solid.size(); i++){
        if(vertex == current_solid[i]){
            number = i;
            break;
        }
    }
    
    return number;
}

vector<Double_t> ConvertToCart(vector<Double_t> vector){
    return { vector[0]*sin(vector[1])*cos(vector[2]),
             vector[0]*sin(vector[1])*sin(vector[2]),
             vector[0]*cos(vector[1])};
}

vector<Double_t> ConvertToSpherical(vector<Double_t> vector){
    TVector3 *vect = new TVector3(vector[0], vector[1], vector[2]);  
 
    return {vect->Mag(), vect->Theta(), vect->Phi()};
}
