#include<iostream>
#include <cmath>
#include <vector>
#include "preparation.hpp"

using namespace std;

class Isingsystem_cubic: public Isingsystem{
private:
      vector<int> system_size;
      void setup_NN(){
          for(int site_idx = 0; site_idx < n_spins; site_idx ++){
              vector<int> r = lattice_coordinate(site_idx);
              spin[site_idx].set_NN(0,site_index(shift_pos_x(r)));
              spin[site_idx].set_NN(1,site_index(shift_pos_y(r)));
              spin[site_idx].set_NN(2,site_index(shift_neg_x(r)));
              spin[site_idx].set_NN(3,site_index(shift_neg_y(r)));
              spin[site_idx].set_NN(4,site_index(shift_pos_z(r)));
              spin[site_idx].set_NN(5,site_index(shift_neg_z(r)));
          }
        
      }

public:
    Isingsystem_cubic(const vector<int> system_size_spec):           // input the three length of the lattice
        Isingsystem(system_size_spec[0] * system_size_spec[1] * system_size_spec[2]),   // calculate the total number of spins
        system_size(system_size_spec){
        set_NN_dim(6);
        setup_NN();
        }
    ~Isingsystem_cubic(){}

    int site_index(const vector<int> lattice_coordinate) const{
        int index = 0;   // find the related index for a coordinate
        index += lattice_coordinate[2] * system_size[0] * system_size[1];
        index += lattice_coordinate[1] * system_size[0] + lattice_coordinate[0];
        return index;
    };

    vector<int> lattice_coordinate(int site_index) const{
        vector<int>coordinate = {0, 0, 0};
        int xy = site_index % (system_size[0] * system_size[1]);
        coordinate[2] = (site_index - xy) / (system_size[0] * system_size[1]);
        coordinate[1] = (xy - xy % system_size[0]) / system_size[0];
        coordinate[0] = xy % system_size[0];
        return coordinate;
    };

    vector<int> shift_pos_x(const vector<int> r_spec) const {                    // 寻找一个坐标对应周边的6个坐标
        vector<int> r(r_spec);
        r[0] = (r[0] + 1) % system_size[0];
        return r;
    };
    vector<int> shift_pos_y(const vector<int> r_spec) const {                    // 寻找一个坐标对应周边的6个坐标
        vector<int> r(r_spec);
        r[1] = (r[1] + 1) % system_size[1];
        return r;
    };
    vector<int> shift_pos_z(const vector<int> r_spec) const {                    // 寻找一个坐标对应周边的6个坐标
        vector<int> r(r_spec);
        r[2] = (r[2] + 1) % system_size[2];
        return r;
    };
    vector<int> shift_neg_x(const vector<int> r_spec) const {                    // 寻找一个坐标对应周边的6个坐标
        vector<int> r(r_spec);
        r[0] = (r[0] - 1 + system_size[0]) % system_size[0];
        return r;
    };
    vector<int> shift_neg_y(const vector<int> r_spec) const {                    // 寻找一个坐标对应周边的6个坐标
        vector<int> r(r_spec);
        r[1] = (r[1] - 1 + system_size[1]) % system_size[1];
        return r;
    };
    vector<int> shift_neg_z(const vector<int> r_spec) const {                    // 寻找一个坐标对应周边的6个坐标
        vector<int> r(r_spec);
        r[2] = (r[2] - 1 + system_size[2]) % system_size[2];
        return r;
    };
};