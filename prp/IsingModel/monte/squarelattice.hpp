#include<iostream>
#include <cmath>
#include <vector>
#include "preparation.hpp"

class Isingsystem_square : public Isingsystem{
private:
    vector<int> system_size;
    void setup_NN(){
        for(int site_idx = 0; site_idx < n_spins; site_idx ++){
            vector<int> r = lattice_coordinate(site_idx);
            spin[site_idx].set_NN(0,site_index(shift_pos_x(r)));
            spin[site_idx].set_NN(1,site_index(shift_pos_y(r)));
            spin[site_idx].set_NN(2,site_index(shift_neg_x(r)));
            spin[site_idx].set_NN(3,site_index(shift_neg_y(r)));
        }
    }

public:
    Isingsystem_square(const vector<int> system_size_spec):
        Isingsystem(system_size_spec[0] * system_size_spec[1]),    // 总晶格数
        system_size(system_size_spec){
        set_NN_dim(4);
        setup_NN();
        }
    ~Isingsystem_square(){}

    int site_index(const vector<int> lattice_coordinate) const{
        return system_size[0] * lattice_coordinate[1] + lattice_coordinate[0];   // 由坐标找到对应的号码，如（3,4）对应27
    };
    vector<int> lattice_coordinate(int site_index) const{
        return {site_index % system_size[0] , (site_index - site_index % system_size[0]) / system_size[0]};  // 由号码找到对应的坐标，如27对应（3,4）
    };

    vector<int> shift_pos_x(const vector<int> r_spec) const {                    // 寻找一个坐标对应周边的四个坐标
        vector<int> r(r_spec);
        r[0] = (r[0] + 1) % system_size[0];
        return r;
    };
    vector<int> shift_pos_y(const vector<int> r_spec) const {                    // 寻找一个坐标对应周边的四个坐标
        vector<int> r(r_spec);
        r[1] = (r[1] + 1) % system_size[1];
        return r;
    };
    vector<int> shift_neg_x(const vector<int> r_spec) const {                    // 寻找一个坐标对应周边的四个坐标
        vector<int> r(r_spec);
        r[0] = (r[0] - 1 + system_size[0]) % system_size[0];
        return r;
    };
    vector<int> shift_neg_y(const vector<int> r_spec) const {                    // 寻找一个坐标对应周边的四个坐标
        vector<int> r(r_spec);
        r[1] = (r[1] - 1 + system_size[1]) % system_size[1];
        return r;
    };
};