#pragma once

namespace online_visualization {
    enum CommandType {
        _null_ = 0,
        _exit_,
        help,
        load,
        list,
        switch_env,
        replay,
        log,
        save,
        batch,
        load_genome,
        set_folder,
        set_file_name,
        set_max_vertices,
        fill_pos,
        clear_pos,
        draw_vertex,
        draw_edge,
        draw_position,
        show_position,
        draw_contig,
        draw_part_of_genome,
        print_paths,
        print_contigs_stats
    };
}

