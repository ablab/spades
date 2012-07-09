#pragma once

#include "all_commands.hpp"
#include "drawing_commands.hpp"
#include "setting_commands.hpp"
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
                //{ "set_folder",         set_folder       }, 
                //{ "set_filename",       set_filename     }, 
                { "set_max_vertices",   set_max_vertices }, 
                { "fill_pos",           fill_pos         }, 
                { "clear_pos",          clear_pos        },
                { "vertex",             draw_vertex      }, 
                { "edge",               draw_edge        }, 
                { "position",           draw_position    } 
                //{ "paths",              paths            }
            };
            
            //VERIFY(info.size() == CommandType::total)
            return CommandNameMapping(info.begin(), info.end());
        }

        static const CommandNameMapping& CommandNameInfo() {
            static CommandNameMapping command_name_mapping = FillCommandNameMapping();
            return command_name_mapping;
        }

        void AddMapping(CommandMapping& mapping, Command* command) {
            mapping.insert(make_pair(command->command_id(), shared_ptr<Command>(command)));
        }

        CommandMapping FillCommandMapping() {
            CommandMapping mapping;
            AddMapping(mapping, new NullCommand);
            AddMapping(mapping, new LoadCommand);
            AddMapping(mapping, new ExitCommand);
            AddMapping(mapping, new ListCommand);
            AddMapping(mapping, new SetMaxVertCommand);
            AddMapping(mapping, new FillPositionCommand);
            AddMapping(mapping, new ListCommand);
            return mapping;
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
            it = command_impl::CommandNameInfo().left.find("null");
        }
        return it->second;
    }

    Command& GetCommand(CommandType command_id) {
        static CommandMapping mapping = command_impl::FillCommandMapping();
        auto it = mapping.find(command_id);
        if (it == mapping.end()) {
            cerr << "No such command, try again" << endl;
            it = mapping.find(CommandType::_null_);
        }
        return *(it->second);
    }

}
