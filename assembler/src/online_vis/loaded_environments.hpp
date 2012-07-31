#pragma once

#include "environment.hpp"

namespace online_visualization {

    typedef shared_ptr<Environment> EnvironmentPtr;

    typedef map<string, EnvironmentPtr> LoadedEnvironments;

    static vector<string>& GetHistory() {
        static vector<string> history;
        return history;
    }

    //static LoadedEnvironments& GetLoadedEnvironments() {
        //static LoadedEnvironments loaded_environments;
        //return loaded_environments;
    //}
}
