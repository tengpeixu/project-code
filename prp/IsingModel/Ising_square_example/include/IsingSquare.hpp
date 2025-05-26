#ifndef isingsystem_square_hpp
#define isingsystem_square_hpp

#include <iostream>
#include <cmath>
#include <cassert>
#include <cstdio>
#include <vector>

#include "IsingSystem.hpp"
using namespace std;

class isingspinonlattice : public IsingSpin
{
private:
    vector<int> position;
    vector<int> NN;

public:
    isingspinonlattice()
    {
        set_dim(1);
        NN = {-1};
    }
    ~isingspinonlattice(){};

    // 函数实现如下
    void set_dim(int dim) { position.assign(dim, 0); }; // 位置向量的大小
    vector<int> _position() const { return position; }
    vector<int> _NN() const { return NN; }
    int _NN(const int bond_idx) const { return NN[bond_idx]; }
    void set_NN(const int bond_idx, const int site_idx) { NN[bond_idx] = site_idx; };
};

class Isingsystem 
{
protected:
    const double J;
    const int n_spins;
    const long long maxrep_state;
    vector<isingspinonlattice> spin;

public:
    Isingsystem(const int n_spins_spec) : 
    J(1.0), 
    n_spins(n_spins_spec),
    maxrep_state(static_cast<long long>(pow(2, n_spins)) - 1)
    {
        spin.resize(n_spins);
    };
    virtual ~Isingsystem(){};

    void set_dim(int dim)
    {
        for (auto &each : spin) //设置每一个（each）的维度
            each.set_dim(dim);
    };

    vector<int> _spin_position(const int site_idx) const { return spin[site_idx]._position(); };
    int _spin_NN(const int site_idx, const int bond_idx) const { return spin[site_idx]._NN(bond_idx); };

    int _n_spins_() const { return n_spins; };

    long long _maxrep_state() const { return maxrep_state; };


    //这一段不再赘述，和一维很像
    int _sz(const int site_idx) const { return spin[site_idx]._sz(); } 
    void set_up_spin(const int site_idx) { spin[site_idx].set_up(); }; 
    void set_dw_spin(const int site_idx) { spin[site_idx].set_dw(); }; 
    void set_spin(const int site_idx, int s_spec)
    {
        spin[site_idx].set_sz(s_spec);
    };
    void flip_spin(const int site_idx) { spin[site_idx].flip(); }; 

    void set_state_by_num(long long rep_state) //和之前一样的
    {
        int i = 0;
        while (rep_state > 0)
        {
            int remainder = rep_state % 2;  
            set_spin(i, 2 * remainder - 1); 
            rep_state /= 2;                 
            i++;
        }

        while (i < n_spins)
        {
            set_spin(i, -1);
            i++;
        }
    }
};

class ising_system_square : public Isingsystem
{
private:
    vector<int> system_size;
    void setup_NN()
    {
        for (int site_idx = 0; site_idx < n_spins; site_idx++)
        {
            vector<int> r = lattice_coordinate(site_idx);
            spin[site_idx].set_NN(0, site_index(shift_pos_x(r)));
            spin[site_idx].set_NN(1, site_index(shift_pos_y(r)));
            spin[site_idx].set_NN(2, site_index(shift_neg_x(r)));
            spin[site_idx].set_NN(3, site_index(shift_neg_y(r)));
        }
    }

public:
    ising_system_square(const vector<int> sys_size) : 
    Isingsystem(calc_n_spins(sys_size)), 
    system_size(sys_size)
    {
        setup_NN();
    };

    int calc_n_spins(const vector<int> &s_size)
    {
        int n_s = 1;
        for (auto &each : s_size)
            n_s *= each;
        return n_s;
    }
    ~ising_system_square(){};

    int site_index(vector<int> lattice_coordinate) const
    {
        int _site_idx = lattice_coordinate[0];
        int w = 1;
        for (int i = 1; i < static_cast<int>(lattice_coordinate.size()); i++)
        {
            for (int j = 0; j < i; j++)
            {
                w *= system_size[j];
            }
            _site_idx = _site_idx + w * lattice_coordinate[i];
            w = 1;
        };
        return _site_idx;
    };

    vector<int> lattice_coordinate(int site_idx) const
    {
        vector<int> l_c(system_size.size(), 0);
        int w = 1;
        for (int i = static_cast<int>(system_size.size() - 1); i >= 0; i--)
        {
            for (int j = 0; j < i; j++)
            {
                w *= system_size[j];
            }
            l_c[i] = site_idx / w;
            site_idx = site_idx % w;
            w = 1;
        };
        return l_c;
    };

    //下面设定周围粒子函数
    vector<int> shift_pos_x(const vector<int> r_spec) const
    {
        vector<int> r(r_spec);
        r[0] = (r[0] + 1) % system_size[0];
        return r;
    }
    vector<int> shift_neg_x(const vector<int> r_spec) const
    {
        vector<int> r(r_spec);
        r[0] = ((r[0] - 1) % system_size[0] < 0) ? ((r[0] - 1) % system_size[0] + system_size[0]) : ((r[0] - 1) % system_size[0]);
        return r;
    }
    vector<int> shift_pos_y(const vector<int> r_spec) const
    {
        vector<int> r(r_spec);
        r[1] = (r[1] + 1) % system_size[1];
        return r;
    }
    vector<int> shift_neg_y(const vector<int> r_spec) const
    {
        vector<int> r(r_spec);
        r[1] = ((r[1] - 1) % system_size[1] < 0) ? ((r[1] - 1) % system_size[1] + system_size[1]) : ((r[1] - 1) % system_size[1]);
        return r;
    }

    int NN(const int site_idx, const int bond_idx) const
    {
        return spin[site_idx]._NN(bond_idx);
    }

    void set_state(vector<bool> state_num)
    {
        for (int site_idx = 0; site_idx < n_spins; site_idx++)
        {
            set_spin(site_idx, 2 * state_num[site_idx] - 1);
        }
    }

    int boundary_length()
    {
        int l = 0;
        for (int i = 0; i < n_spins; i++)
        {
            vector<int> a = lattice_coordinate(i);
            int N0 = site_index(shift_pos_x(a));
            int N1 = site_index(shift_pos_y(a));
            int N2 = site_index(shift_neg_x(a));
            int N3 = site_index(shift_neg_y(a));
            l = l - (_sz(i) * _sz(N0) - 1) / 2 - (_sz(i) * _sz(N1) - 1) / 2 - (_sz(i) * _sz(N2) - 1) / 2 - (_sz(i) * _sz(N3) - 1) / 2;
        }
        l = l / 2;
        return l;
    }

    double eval_energy()
    {
        double e = J * (boundary_length() - (2 * n_spins - boundary_length()));
        return e;
    }

    double eval_Mz()
    {
        double mz = 0;
        for (int i = 0; i < n_spins; i++)
        {
            mz = mz + _sz(i);
        }
        return mz;
    }
};

#endif