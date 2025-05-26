#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <iostream>
#include <vector>
#include "IsingSquare.hpp"
#include "IsingSystem.hpp"

TEST_CASE("Ising_square", "coordinate->site_idx") 
{
    const vector <int> sys_size={6,6};
    ising_system_square model(sys_size);

    SECTION("coordinate->site_idx") 
    {
        const vector<int> l_c= {3,4};
        const vector<int> l_cx= {3,3};
        REQUIRE(model.site_index(l_c)==27);
        REQUIRE(model.lattice_coordinate(27)== l_c);
    };
    SECTION("NN")
    {
        constexpr int i = 21;
        REQUIRE(model.NN(i,0)==22);
        REQUIRE(model.NN(i,1)==27);
        REQUIRE(model.NN(i,2)==20);
        REQUIRE(model.NN(i,3)==15);
    }

    SECTION("pi state")
    {
    vector <bool> state({1,1,0,1,1,1,0,0,1,1,1,0,1,1,1,1,0,1,0,0,0,0,0,0,1,1,0,1,0,1,1,1,0,0,0,0});
    model.set_state(state);
    REQUIRE(model.eval_Mz()==2);
    REQUIRE_THAT(model.eval_energy(),Catch::Matchers::WithinULP(-4.0,4));
    }
};

TEST_CASE("Ising_square 3D", "coordinate->site_idx,3D") 
{
    const vector <int> sys_size1={3,3,3};
    ising_system_square model1(sys_size1);

    SECTION("coordinate->site_idx") 
    {
        const vector<int> l_c1= {1,1,1};
        const vector<int> l_cx1= {1,1,2};
        REQUIRE(model1.site_index(l_c1)==13);
        REQUIRE(model1.lattice_coordinate(13)== l_c1);
    };
};

TEST_CASE("shift", "") 
{
    const vector <int> sys_size2={6,6};
    ising_system_square model2(sys_size2);

    SECTION("coordinate->site_idx") 
    {
        const vector<int> l_c2= {3,3};
        const vector<int> l_c21= {0,3};
        const vector<int> l_c2_n_x= {2,3};
        const vector<int> l_c2_p_x= {4,3};
        const vector<int> l_c2_n_y= {3,2};
        const vector<int> l_c2_p_y= {3,4};
        const vector<int> l_c21_n_x= {5,3};
        REQUIRE(model2.shift_neg_x(l_c2)==l_c2_n_x);
        REQUIRE(model2.shift_pos_x(l_c2)==l_c2_p_x);
        REQUIRE(model2.shift_neg_y(l_c2)==l_c2_n_y);
        REQUIRE(model2.shift_pos_y(l_c2)==l_c2_p_y);
        REQUIRE(model2.shift_neg_x(l_c21)==l_c21_n_x);
    };
};