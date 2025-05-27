#include<iostream>
#include <cassert>
#include <cmath>
#include <bitset>
#include <string>
#include <vector>
# pragma once

using namespace std;

class Isingspin{
private:
    int sz; // to decide the ising spin situation

public:
    int _sz() const { return sz;} // const means no change to class attributes
    void set_up(){sz = 1;};
    void set_down(){sz = -1;};
    void set_sz(int sz_spec){
        assert(sz_spec == 1 || sz_spec == -1);
        sz = sz_spec;
    };
    void flip(){sz *= -1;};
};

class IsingSpinOnLattice : public Isingspin{
private:
    vector<int> position;
    vector<int> NN;

public:
    IsingSpinOnLattice(){
        NN = {-1};        
    };
    ~IsingSpinOnLattice(){}

    vector<int> _position() const{return position;};
    vector<int> _NN() const{return NN;};
    int _NN(const int bond_idx) const{return NN[bond_idx];};
    void set_NN(const int bond_idx, const int site_idx){
        NN[bond_idx] = site_idx;
    }
    void set_NN_dim(int dim){NN.resize(dim);};
};

class Observable {
private:
    double z;
    double q;
    double q_2;
    double q_4;

public:
    Observable(){reset();};
    ~Observable(){};

    void reset(){z = q = q_2 = q_4 = 0;}
    
    void update_direct(const double& data_sample, const double data_W){
        z += data_W;
        q += data_W * data_sample;
        q_2 += data_W * data_sample * data_sample;
        q_4 += data_sample * data_sample * data_sample * data_sample * data_W;
    }

    void normalize_direct(){
        q = q / z;
        q_2 /= z;
        q_4 /= z;
    }
        
    double _z() const {return z;};
    double _q() const {return q;};
    double _q2() const {return q_2;};
    double _q4() const {return q_4;};
    double _q_fluc() const {return q_2 - q * q;};
    double _q_binder() const {return q_4 / (q_2 * q_2);};
};


class Isingsystem{
protected:
    const double J;
    const int n_spins;
    const long long maxrep_state;
    int NN_dim = 1;
    vector<IsingSpinOnLattice> spin;

public:
    Isingsystem(const int n_spins_spec): J(1.0), n_spins(n_spins_spec),                    // input the total number of spins
        maxrep_state(static_cast<long long>(pow(2, n_spins))-1){
            spin.resize(n_spins);
        };
    virtual ~Isingsystem(){};
    
    void set_NN_dim(int dim){
        for (auto& each:spin) each.set_NN_dim(dim);
        NN_dim = dim;
        };                       // set the dimension for all the spins(3)
    int _n_spins(){return n_spins;};
    float _J(){return J;};
    vector<int> _spin_position(const int site_idx) const{return spin[site_idx]._position();};          // return one position(position and coordinate are different)
    vector<int> _spin_NN(const int site_idx) const{return spin[site_idx]._NN();};
    int _spin_NN(const int site_idx, const int bond_idx) const {return spin[site_idx]._NN(bond_idx);};

    void set_state(long long rep_state){
        bitset<64> binary(rep_state);
        string binary_str = binary.to_string();
        string padded_binary_str = binary_str.substr(binary_str.length() - n_spins);
        for(int i = 0; i < n_spins; i++){
            if(padded_binary_str[i] == '1'){spin[i].set_up();}
            if(padded_binary_str[i] == '0'){spin[i].set_down();}
        }
    };
    
    int eval_mz() const{                     // 计算M
        int mz = 0;
        for(int i = 0; i < n_spins; i++){
            mz += spin[i]._sz();
        }
        return mz;
    };

    int NN(const int site_idx, const int bond_idx) const {                
        return spin[site_idx]._NN(bond_idx);
    };

    double eval_h() const{                      // 计算H
        float h = 0;
        for(int i = 0; i < n_spins; i++){
            for(int j = 0; j < NN_dim; j++){
                h = h + (-0.5 * J * spin[i]._sz() * spin[NN(i, j)]._sz()); 
            };
        };
        return h;
    };

    vector<vector<double>> evaluation(vector<double> temp){           

        vector<double> M2(temp.size(), 0);
        vector<double> H(temp.size(), 0);
        vector<double> H2(temp.size(), 0);
        vector<double> M4(temp.size(), 0);  
        vector<double> C(temp.size(), 0);
        vector<double> binder_ratio(temp.size(), 0);
        vector<Observable> h(temp.size());
        vector<Observable> mz(temp.size());                      
        
        for(long long state = 0; state < maxrep_state + 1; state ++){
            set_state(state);
            double hh = eval_h();
            double mag = eval_mz();
            for(int i = 0; i < temp.size(); i++){
                h[i].update_direct(hh, exp(- hh/temp[i]));
                mz[i].update_direct(mag, exp(- hh/temp[i]));
            }
        }
        
        for(int i = 0; i < temp.size(); i++){
            h[i].normalize_direct();
            mz[i].normalize_direct();
            M2[i] = mz[i]._q2();
            H[i] = h[i]._q();
            C[i] = h[i]._q_fluc();
            binder_ratio[i] = mz[i]._q_binder();
        }
        
        return {M2, H, C, binder_ratio};
    };
};
