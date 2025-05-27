#include<iostream>
#include"cubiclattice.hpp"
#include"kagomelattice.hpp"
#include"squarelattice.hpp"

using namespace std;

int main()
{   vector<double> temp(120); 
    double start = 0.05;
    for (int j = 0; j < 120; j++){
        temp[j] = start;
        start += 0.05;
    }
    
    // cubic lattice
    int i1 = 2;
    int i2 = 2;
    int i3 = 2;
    Isingsystem_cubic model ({i1, i2, i3});

    vector<vector<double>> result_vector = model.evaluation(temp);
    for(int j = 0; j < temp.size(); j++){
        double C = (result_vector[2][j]) / (temp[j] * temp[j] * i1 * i2 * i3);
        double M2 = result_vector[0][j] / (i1 * i2 * i3 * i1 * i2 * i3);
        double E = result_vector[1][j] / (i1 * i2 * i3);
        double binder_ratio = result_vector[3][j];
        cout << "2,2,2" << "\t" << temp[j] << "\t" << C << "\t" << M2 << "\t" << E << "\t" << binder_ratio << endl;
    }

    // kagome lattice
    for(int i = 2; i<11;){
        Isingsystem_kagome model ({i,i});
        
        string sstr = "./Ising_kagome_" + std::to_string(i) + ".txt";
        ofstream output;

        output.open(sstr, ios::out | ios::trunc);
        result_vector = model.evaluation(temp);
        for(int j = 0; j < size; j++){
            int number = (i * i - (i / 2 ) * ( i / 2 ));
            double C = (result_vector[2][j] - result_vector[1][j] * result_vector[1][j]) / (temp[j] * temp[j] * number);
            double M2 = result_vector[0][j] / (pow(number , 2));
            double E = result_vector[1][j] / number;
            double U = result_vector[3][j];
            output << i << "\t" << temp[j] << "\t" << C << "\t" << M2 << "\t" << E << "\t" << U << endl;            
        
        }
        output.close();
        cout << i << "*" << i << "已完成！" << endl;
        i += 2;
     }

    // square lattice
    for(int i = 2; i< 6; i++){
        Isingsystem_square model ({i, i});
        result_vector = model.evaluation(temp);
        for(int j = 0; j < temp.size(); j++){
            double C = (result_vector[2][j]) / (temp[j] * temp[j] * i * i);
            double M2 = result_vector[0][j] / pow(i, 4);
            double E = result_vector[1][j] / pow(i, 2);
            double binder_ratio = result_vector[3][j];
            cout << "2,2,2" << "\t" << temp[j] << "\t" << C << "\t" << M2 << "\t" << E << "\t" << binder_ratio << endl;
        }
    }
}