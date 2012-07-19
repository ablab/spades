#pragma once

namespace online_visualization {
    enum CommandType {
        _null_ = 0,
        _exit_,
        load,
        list,
        switch_env,
        replay,
        set_folder,
        set_file_name,
        set_max_vertices,
        fill_pos,
        clear_pos,
        draw_vertex,
        draw_edge,
        draw_position,
        print_paths
        //total = 2
    };
}

