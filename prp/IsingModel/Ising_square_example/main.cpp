#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include "IsingSquare.hpp"
#include "IsingSystem.hpp"
using namespace std;

int main()
{   
    int a;
    cout << "请输入正方形边长：" ;
    cin >> a;
    cin.get();
    string sstr = "./ising_square_dataL" + std::to_string(a) + ".txt";
    ofstream fout(sstr, ios::trunc);
    fout << "T/|J|"
         << " "
         << "M2/N2"
         << " "
         << "H/N"
         << " "
         << "C" 
         << endl;

    vector<int> sys_size = {a, a};
    ising_system_square model(sys_size);
    model.set_state_by_num(0);
    double N = model._n_spins_();
    double H0 = model.eval_energy();
    long long max_state = model._maxrep_state();

    vector<double> energy_values(max_state + 1, 0.0);
    vector<double> Mz_values(max_state + 1, 0.0);

    for (int i = 0; i <= max_state; i++)
    {
        model.set_state_by_num(i);
        energy_values[i] = model.eval_energy();
        Mz_values[i] = model.eval_Mz();
    }

    for (int j = 1; j <= 80; j++)
    {
        double T = 0.05 * j;
        double z = 0, H = 0, H2 = 0, C = 0, M2 = 0;

        for (int i = 0; i <= max_state; i++)
        {
            z += exp(-(energy_values[i] - H0) / T);
            H += energy_values[i] * exp(-(energy_values[i] - H0) / T);
            H2 += pow(energy_values[i], 2) * exp(-(energy_values[i] - H0) / T);
            M2 += pow(Mz_values[i], 2) * exp(-(energy_values[i] - H0) / T);
        }

        H = H / z;
        H2 = H2 / z;
        C = (H2 - pow(H, 2)) / pow(T, 2);
        M2 = M2 / z;
        fout << T << " " << M2 / pow(N, 2) << " " << H / N << " " << C << endl;
    }

    fout.close();
    cout << "目标已完成！" << endl;
    return 0;
}