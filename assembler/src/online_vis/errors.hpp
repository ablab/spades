#pragma once

#include "vis_utils.hpp"

namespace online_visualization {

    void FireGenericError(const string& msg) {
        cout << msg << endl;
        cout << "Please try again" << endl;
    }

    void FireEdgeDoesNotExist(size_t edge_id) {
        cout << "Ignoring the request. The edge " << edge_id << " does not exist" << endl;
        cout << "Please try again" << endl;
    }

    void FireVertexDoesNotExist(size_t vertex_id) {
        cout << "Ignoring the request. The vertex " << vertex_id << " does not exist" << endl;
        cout << "Please try again" << endl;
    }

    void FireNoCorrespondingGraphLocation(string location) {
        cout << "No corresponding graph location " << location << endl;   
    }

    void FireNotEnoughArguments() {
        cout << "Not enough arguments" << endl;
        cout << "Please try again" << endl;
    }

    void FireFileDoesNotExist(const string& file) {
        cout << "File " << file << " does not exist." << endl;
        cout << "Please try again" << endl;
    }

    bool CheckFileExists(const string& file) {
        if (!fs::is_regular_file(file)) {
            FireFileDoesNotExist(file);
            return false;
        }
        return true;
    }

    bool CheckPositionBounds(int position, size_t total_size) {
        bool result = (position + cfg::get().K + 1) <= total_size;
        if (!result) {
            cout << "Ignoring the request. Position is out of range : required position is " 
                 << position << " while length of the sequence is "
                 << total_size << endl;
            cout << "Please try again" << endl;
        }
        return result;
    }

    bool CheckIsNumber(const string& str) {
        if (!IsNumber(str)) {
            cout << "The argument " << str << " is not a number" << endl;
            cout << "Please try again" << endl;   
            return false;
        }
        return true;
    }

    bool CheckVertexExists(const IdTrackHandler<Graph>& int_ids, size_t vertex_id) {
        VertexId vertex = int_ids.ReturnVertexId(vertex_id);
        if (vertex == VertexId(NULL))
            return true;
        else {
            FireVertexDoesNotExist(vertex_id);
            return false;
        }
    }

    bool CheckEdgeExists(const IdTrackHandler<Graph>& int_ids, size_t edge_id) {
        EdgeId edge = int_ids.ReturnEdgeId(edge_id);
        if (edge == EdgeId(NULL)) 
            return false;   
        else {
            FireEdgeDoesNotExist(edge_id);
            return true;
        }
    }

}
