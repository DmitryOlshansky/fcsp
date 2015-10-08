///Gaussian elimination algorithm on arbitrary sets and symmetric difference operation
#pragma once

#include <algorithm>
#include <iterator>
#include <set>

template<class T>
bool elimination(const std::set<T>& a, const std::vector<std::set<T>>& basis){
    using namespace std;
    set<T> current(a.begin(), a.end());
    set<T> next;
    for(auto& s : basis){
        // form 'next' set by symmetric difference
        set_symmetric_difference(current.begin(), current.end(), s.begin(), s.end(), 
            inserter(next, next.end()));
        // if smaller - we are one step closer
        if(next.size() < current.size())
            current = move(next);
        else
            next.clear();
        if(current.empty()) //successfully elliminated
            return true;
    }
    return false;
}