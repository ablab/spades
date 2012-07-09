#pragma once

namespace online_visualization {
    enum CommandType {
        _exit_ = 0,
        load
        //set_folder,
        //set_filename,
        //set_max_vertices,
        //fill_pos,
        //clear_pos,
        //vertex,
        //edge,
        //position,
        //paths,
        //total = 2
    };

    namespace command_impl {

        typedef online_visualization::CommandType CommandType;

        typedef bimap<string, CommandType> CommandNameMapping;

        static const CommandNameMapping FillCommandNameMapping() {
            vector<CommandNameMapping::value_type> info = 
            { 
                { "exit",               _exit_           }, 
                { "load",               load             }
                //{ "set_folder",         set_folder       }, 
                //{ "set_filename",       set_filename     }, 
                //{ "set_max_vertices",   set_max_vertices }, 
                //{ "fill_pos",           fill_pos         }, 
                //{ "clear_pos",          clear_pos        }, 
                //{ "vertex",             vertex           }, 
                //{ "edge",               edge             }, 
                //{ "position",           position         }, 
                //{ "paths",              paths            }
            };
            
            //VERIFY(info.size() == CommandType::total)
            return CommandNameMapping(info.begin(), info.end());
        }

        static const CommandNameMapping& CommandNameInfo() {
            static CommandNameMapping command_name_mapping = FillCommandNameMapping();
            return command_name_mapping;
        }


    }


    //static const string& CommandName(CommandType command) {

        //auto it = command_impl::CommandNameInfo().right.find(command);

        //VERIFY_MSG(it != command_impl::CommandNameInfo().right.end(),
                //"No name for command id = " << command);
        //return it->second;
    //}

    static CommandType CommandId(string name) {
        auto it = command_impl::CommandNameInfo().left.find(name);
        VERIFY_MSG(it != command_impl::CommandNameInfo().left.end(),
                "There is no command type with name = " << name);

        return it->second;
    }

}
