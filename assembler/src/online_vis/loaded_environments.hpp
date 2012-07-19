#pragma once

#include "environment.hpp"

namespace online_visualization {

    typedef map<string, shared_ptr<Environment> > LoadedEnvironments;


    static vector<string>& GetHistory() {
        static vector<string> history;
        return history;
    }

    static LoadedEnvironments& GetLoadedEnvironments() {
        static LoadedEnvironments loaded_environments;
        return loaded_environments;
    }
}
