#pragma once

#include "command_type.hpp"

namespace online_visualization {

    typedef map<CommandType, shared_ptr<Command> > CommandMapping;

    namespace command_impl {

        typedef bimap<string, CommandType> CommandNameMapping;

        static const CommandNameMapping FillCommandNameMapping() {
            vector<CommandNameMapping::value_type> info = 
            { 
                { "null",               _null_           }, 
                { "exit",               _exit_           }, 
                { "load",               load             },
                { "list",               list             },
                { "switch",             switch_env       },
                { "rep",                replay           },
                { "set_folder",         set_folder       }, 
                { "set_file_name",      set_file_name    }, 
                { "set_max_vertices",   set_max_vertices }, 
                { "fill_pos",           fill_pos         }, 
                { "clear_pos",          clear_pos        },
                { "vertex",             draw_vertex      }, 
                { "edge",               draw_edge        }, 
                { "position",           draw_position    }, 
                { "paths",              print_paths      }
            };
            
            //VERIFY(info.size() == CommandType::total)
            return CommandNameMapping(info.begin(), info.end());
        }

        static const CommandNameMapping& CommandNameInfo() {
            static CommandNameMapping command_name_mapping = FillCommandNameMapping();
            return command_name_mapping;
        }

        void AddMapping(CommandMapping& mapping, shared_ptr<Command>& command) {
            mapping.insert(make_pair(command->command_id(), command));
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
        if (it == command_impl::CommandNameInfo().left.end()) {
            cout << "No such command, try again" << endl;
            it = command_impl::CommandNameInfo().left.find("null");
        }
        return it->second;
    }

    static CommandMapping& GetCommandMapping() {
        static CommandMapping mapping;
        return mapping;
    }


    Command& GetCommand(CommandType command_id) {
        const CommandMapping& mapping = GetCommandMapping();
        auto it = mapping.find(command_id);
        return *(it->second);
    }

    void AddCommand(shared_ptr<Command> command) {
        CommandMapping& mapping = GetCommandMapping();
        command_impl::AddMapping(mapping, command);
    }

}
