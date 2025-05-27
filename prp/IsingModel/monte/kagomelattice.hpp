#include <iostream>
#include <cmath>
#include <algorithm>
#include "preparation.hpp"

using namespace std;


class Isingsystem_kagome : public Isingsystem
{
private:
    vector<int> system_size;
    void setup_NN()
    {
        for (int i = 0; i < n_spins; i++)
        {
            vector<int> r = lattice_coordinate(i);
            spin[i].set_NN(0, surr(r)[0]);
            spin[i].set_NN(1, surr(r)[1]);
            spin[i].set_NN(2, surr(r)[2]);
            spin[i].set_NN(3, surr(r)[3]);
        }
    }

public:
    Isingsystem_kagome(const vector<int> system_size_spec) : Isingsystem((system_size_spec[0] * system_size_spec[1]-(system_size_spec[1]/2)* (system_size_spec[0]/2))), // 总晶格数
                                                             system_size(system_size_spec)
    {
        setup_NN();
        set_NN_dim(4);
    }
    ~Isingsystem_kagome() {}

    int _total_number(){return system_size[1] * system_size[0] - (system_size[0]/2) * (system_size[1]/2);};

    vector<int> _system_size(){return system_size;};

    int site_index(const vector<int> lattice_coordinate) const
    {
        return (system_size[0] * lattice_coordinate[1] + lattice_coordinate[0] / (1 + lattice_coordinate[1] % 2) - (system_size[0] / 2) * (lattice_coordinate[1] / 2)); // 由坐标找到对应的号码，如（3,4）对应27
    };



    vector<int> lattice_coordinate(int site_index) const
    {
        int b, x, y;
        int tmp = system_size[0] + system_size[0]/2;
        site_index = 1 + site_index;
        b = site_index % tmp;
        if (b == 0)
        {
            x = (site_index / tmp) * 2;
            y = system_size[0]-1;
        }
        else
        {
            if (b % system_size[0] == 0)
            {
                y = system_size[0];
                x = (site_index / tmp) * 2 + 1;
            }
            else
            {
                if (b % system_size[0] == 0)
                {
                    x = (site_index / tmp) * 2 +1;
                    y = system_size[0];
                }
                else
                {
                    x = (site_index / tmp) * 2 +b/system_size[0]+1;
                    y = (b % system_size[0] - 1) * (system_size[1] % 2 + 1) + 1;
                }
            }
        }
        return {y - 1, x - 1};
    };

    vector<int> shift_pos_x(const vector<int> r_spec) const
    { // 寻找一个坐标对应周边的四个坐标
        vector<int> r(r_spec);
        r[0] = (r[0] + 1) % system_size[0];
        return r;
    };
    vector<int> shift_neg_x(const vector<int> r_spec) const
    { // 寻找一个坐标对应周边的四个坐标
        vector<int> r(r_spec);
        r[0] = (r[0] - 1 + system_size[0]) % system_size[0];
        return r;
    };
    vector<int> shift_pos_y(const vector<int> r_spec) const
    { // 寻找一个坐标对应周边的四个坐标
        vector<int> r(r_spec);
        if (r[1] + 1 == system_size[1])
        {
            r[0] = ((system_size[1]) / 2 + r[0])%system_size[0];
            r[1] = 0;
        }
        else
        {
            r[1] = (r[1] + 1) % system_size[1];
        }
        return r;
    };
    vector<int> shift_neg_y(const vector<int> r_spec) const
    { // 寻找一个坐标对应周边的四个坐标
        vector<int> r(r_spec);
        if (r[1] - 1 == -1)
        {
            r[0] = (r[0] - (system_size[1]) / 2 + system_size[0]) % system_size[0];
            r[1] = system_size[1] - 1;
        }
        else
        {
            r[1] = (r[1] - 1 + system_size[1]) % system_size[1];
        }
        return r;
    };
    vector<int> specical_1(const vector<int> r_spec) const
    { // 寻找一个坐标对应周边的四个坐标
        vector<int> r(r_spec);
        if (r[1] == 0)
        {
            r[0] = (r[0] - (system_size[1]) / 2 + 1 + system_size[0]) % system_size[0];
            r[1] = system_size[1] - 1;
        }
        else
        {
            r[0] = (r[0] + 1) % system_size[0];
            r[1] = (r[1] - 1 + system_size[1]) % system_size[1];
        }
        return r;
    };
    vector<int> specical_2(const vector<int> r_spec) const
    { // 寻找一个坐标对应周边的四个坐标
        vector<int> r(r_spec);
        if (r[1] + 1 == system_size[1])
        {
            r[0] = (r[0] + (system_size[1]) / 2 - 1 + system_size[0]) % system_size[0];
            r[1] = 0;
        }
        else
        {
            r[1] = (r[1] + 1) % system_size[1];
            r[0] = (r[0] - 1 + system_size[0]) % system_size[0];
        }
        return r;
    };

    vector<int> surr(const vector<int> r_spec) const
    {
        vector<int> r(r_spec);
        vector<int> tmp{0, 0, 0, 0};
        if (r[0] % 2 == 0)
        {
            if (r[1] % 2 == 1) // x+y+x-y-
            {
                tmp[0] = site_index(specical_1(r));
                tmp[1] = site_index(shift_pos_y(r));
                tmp[2] = site_index(specical_2(r));
                tmp[3] = site_index(shift_neg_y(r));
            }
            else
            {
                tmp[0] = site_index(shift_pos_x(r));
                tmp[1] = site_index(shift_pos_y(r));
                tmp[2] = site_index(shift_neg_x(r));
                tmp[3] = site_index(shift_neg_y(r));
            }
        }
        else
        {
            tmp[0] = site_index(shift_pos_x(r));
            tmp[1] = site_index(specical_2(r));
            tmp[2] = site_index(shift_neg_x(r));
            tmp[3] = site_index(specical_1(r));
        }
        return tmp;
    }
};
/**
 *                                         ,s555SB@@&                          
 *                                      :9H####@@@@@Xi                        
 *                                     1@@@@@@@@@@@@@@8                       
 *                                   ,8@@@@@@@@@B@@@@@@8                      
 *                                  :B@@@@X3hi8Bs;B@@@@@Ah,                   
 *             ,8i                  r@@@B:     1S ,M@@@@@@#8;                 
 *            1AB35.i:               X@@8 .   SGhr ,A@@@@@@@@S                
 *            1@h31MX8                18Hhh3i .i3r ,A@@@@@@@@@5               
 *            ;@&i,58r5                 rGSS:     :B@@@@@@@@@@A               
 *             1#i  . 9i                 hX.  .: .5@@@@@@@@@@@1               
 *              sG1,  ,G53s.              9#Xi;hS5 3B@@@@@@@B1                
 *               .h8h.,A@@@MXSs,           #@H1:    3ssSSX@1                  
 *               s ,@@@@@@@@@@@@Xhi,       r#@@X1s9M8    .GA981               
 *               ,. rS8H#@@@@@@@@@@#HG51;.  .h31i;9@r    .8@@@@BS;i;          
 *                .19AXXXAB@@@@@@@@@@@@@@#MHXG893hrX#XGGXM@@@@@@@@@@MS        
 *                s@@MM@@@hsX#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@&,      
 *              :GB@#3G@@Brs ,1GM@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@B,     
 *            .hM@@@#@@#MX 51  r;iSGAM@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@8     
 *          :3B@@@@@@@@@@@&9@h :Gs   .;sSXH@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@:    
 *      s&HA#@@@@@@@@@@@@@@M89A;.8S.       ,r3@@@@@@@@@@@@@@@@@@@@@@@@@@@r    
 *   ,13B@@@@@@@@@@@@@@@@@@@5 5B3 ;.         ;@@@@@@@@@@@@@@@@@@@@@@@@@@@i    
 *  5#@@#&@@@@@@@@@@@@@@@@@@9  .39:          ;@@@@@@@@@@@@@@@@@@@@@@@@@@@;    
 *  9@@@X:MM@@@@@@@@@@@@@@@#;    ;31.         H@@@@@@@@@@@@@@@@@@@@@@@@@@:    
 *   SH#@B9.rM@@@@@@@@@@@@@B       :.         3@@@@@@@@@@@@@@@@@@@@@@@@@@5    
 *     ,:.   9@@@@@@@@@@@#HB5                 .M@@@@@@@@@@@@@@@@@@@@@@@@@B    
 *           ,ssirhSM@&1;i19911i,.             s@@@@@@@@@@@@@@@@@@@@@@@@@@S   
 *              ,,,rHAri1h1rh&@#353Sh:          8@@@@@@@@@@@@@@@@@@@@@@@@@#:  
 *            .A3hH@#5S553&@@#h   i:i9S          #@@@@@@@@@@@@@@@@@@@@@@@@@A. 
 *
 *
 *    又看源码，看你妹妹呀！
 */

