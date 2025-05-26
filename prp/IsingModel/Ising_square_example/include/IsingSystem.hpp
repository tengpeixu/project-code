#pragma once
#ifndef IsingSystem_hpp
#define IsingSystem_hpp

#include <iostream> 
#include <cassert> 
#include <cmath>
#include <vector>

// 下面我们定义IsingSpin类，表示粒子的状态
class IsingSpin
{
private:
    int sz; /* +/-1 */

public:
    IsingSpin() { sz = 1; };
    ~IsingSpin(){};

    int _sz() const { return sz; }; // 定义粒子状态函数
    void set_up() { sz = 1; };
    void set_dw() { sz = -1; };
    void set_sz(int sz_spec)
    {
        assert(sz_spec == 1 || sz_spec == -1); // 保证sz_spec一定是+-1 //类似自检程序
        sz = sz_spec;
    };
    void flip() { sz *= -1; }; // 定义函数flip（改变正负）
};


#endif