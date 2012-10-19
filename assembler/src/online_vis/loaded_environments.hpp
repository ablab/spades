#pragma once

#include "environment.hpp"

namespace online_visualization {

    template <class Env>
    class LoadedEnvironments : public map<string, shared_ptr<Env> > {};
//    typedef map<string, shared_ptr<Env> > LoadedEnvironments;

    static vector<string>& GetHistory() {
        static vector<string> history;
        return history;
    }

    //static LoadedEnvironments& GetLoadedEnvironments() {
        //static LoadedEnvironments loaded_environments;
        //return loaded_environments;
    //}
}
