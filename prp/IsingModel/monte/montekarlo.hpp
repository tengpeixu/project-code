#include <iostream>
#include "cubiclattice.hpp"
#include <vector>
#include <cmath>
#include <random>
#include <chrono>

class mc_cubic: public Isingsystem_cubic{
private:
    unsigned random_seed = 0;
public:
    mc_cubic(const vector<int> system_size):
        Isingsystem_cubic(system_size){
        }
    ~mc_cubic(){}

    void set_seed(unsigned seed){random_seed = seed;}

    unsigned _random_seed(){return random_seed;}

    void randomize_spins(){                 // initialization
        unsigned seed = random_seed;
        static std::mt19937 gen(seed);
        static std::uniform_real_distribution<> dis(0, 1);
        for(int i=0; i<n_spins; i++){
            bool random_number = dis(gen) < 0.5;
            int spin_state = 2 * random_number - 1;
            spin[i].set_sz(spin_state);
        }
    }

    void ordered_spins(){
        for(int i=0; i<n_spins; i++){
            spin[i].set_up();
        }
    }

    double eval_local_energy(int index){
        double local_energy = 0;
        for(int i=0; i<6; i++){
            local_energy += - (spin[index]._sz() * spin[NN(index, i)]._sz());
        }
        return local_energy;
    }
};


class monte_carlo{
private:
    double temperature;
    mc_cubic state_before;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis;

public:
    monte_carlo(const vector<int> system_size):
        gen(0),
        dis(0,1),
        state_before(system_size){
        state_before.ordered_spins();
    }
    ~monte_carlo(){}

    void set_temperature(double temp){temperature = temp;}

    void randomize_spins(){state_before.randomize_spins();}

    void ordered_spins(){state_before.ordered_spins();}

    unsigned random_seed(){return state_before._random_seed();}

    void set_seed(int randomseed){
        state_before.set_seed(randomseed);
        gen.seed(randomseed);
    }

    void run_mc(){
        static int index = 0;
        index += 1;
        if(index >= state_before._n_spins()){index = 0;}
        state_before.flip_spin(index);
        double Ef = state_before.eval_local_energy(index);
        state_before.flip_spin(index);
        double Ei = state_before.eval_local_energy(index);
        if(Ei >= Ef){
            state_before.flip_spin(index);
        }
        if(Ei < Ef){
            double p = exp(-(Ef - Ei) / temperature);
            double random_number = dis(gen);
            if(random_number < p){state_before.flip_spin(index);}
        }  
    }

    vector<vector<double>> sweep(int sweep_time){
        vector<double> E(sweep_time);
        vector<double> M(sweep_time);
        vector<double> M_2(sweep_time);
        for(int time=0; time<sweep_time; time++){
            run_mc();
            double e = state_before.eval_h();
            double m = state_before.eval_mz();
            E[time] = e;
            M[time] = m;
            M_2[time] = m * m;
        }
        return {E, M, M_2};
    };

    vector<double> statistic(vector<double> X){
        double sum = accumulate(begin(X), end(X), 0.0);  
        double mean =  sum / X.size();
        double accum  = 0.0;  
        int size_X = X.size();
        for (int i = 0; i < size_X; i++)
        {
            accum += (X[i]-mean)*(X[i]-mean);
        }
        double stdev = sqrt(accum/(size_X-1));//标准误差
        return {mean, stdev};
    }

    vector<vector<double>> bins_analysis(int bins, int bins_scale){
        vector<double> E(bins);
        vector<double> M(bins);
        vector<double> M_2(bins);
        sweep(500);
        for(int bin=0; bin<bins; bin++){
            vector<vector<double>> result = sweep(bins_scale);
            double sum = accumulate(begin(result[0]), end(result[0]), 0.0);
            E[bin] = sum / bins_scale / state_before._n_spins();
            double sum1 = accumulate(begin(result[1]), end(result[1]), 0.0);
            M[bin] = sum1 / bins_scale / state_before._n_spins();
            double sum2 = accumulate(begin(result[2]), end(result[2]), 0.0);
            M_2[bin] = sum2 / bins_scale / pow(state_before._n_spins(),2);
            sweep(bins_scale / 3);
        }
        return{E, M, M_2};
    }

};