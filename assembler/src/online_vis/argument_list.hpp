#pragma once

namespace online_visualization {
    
    class ArgumentList {
        private:
            map<string, string> options;
            set<string> short_options;
            vector<string> arguments;
            //vector<string> optional_arguments;

            const vector<string> SplitInTokens(stringstream& args) const { 
                vector<string> answer;
                while (!args.eof()) {
                    string arg;
                    args >> arg;
                    answer.push_back(arg);
                }
                return answer;
            }

            pair<string, string> ParseOption(const string& arg) const {
                string opt_name;
                string opt_value;

                size_t i = 2;
                for (; i < arg.size() && arg[i] != '='; ++i) {
                    opt_name = opt_name + arg[i];   
                }
                for (; i < arg.size(); ++i) {
                    opt_value = opt_value + arg[i];   
                }
            
                if (opt_value == "")
                    opt_value = "true";
                

                return make_pair(opt_name, opt_value);
            }

            vector<string> ParseShortOption(const string& arg) const {
                vector<string> result;
                size_t i = 1;
                for (; i < arg.size(); ++i) {
                    result.push_back((const char*) arg[i]);
                }       
                return result;
            }
        
            void ParseArguments(const vector<string>& args) {
                for (size_t i = 0; i < args.size(); ++i) {
                    if (args[i][0] == '-' && args[i][1] == '-') {
                        //--smth=<smth>
                        pair<string, string> opt_val = ParseOption(args[i]);

                    }
                    else if (args[i][0] == '-') {
                        const vector<string>& short_opt = ParseShortOption(args[i]);
                        short_options.insert(short_opt.begin(), short_opt.end());
                    }
                    else {
                        arguments.push_back(args[i]);    
                    }
                }
            }

        public:
            ArgumentList() {
            }

            ArgumentList(stringstream& stream) {
                const vector<string>& args = SplitInTokens(stream);
                ParseArguments(args);
            }

            string operator[](const string& option_name) const {
                // usual option
                if (options.find(option_name) == options.end()) {
                    return "null";   
                }
                string result = options.find(option_name)->second;
                return result;
            }

            bool contains(string short_opt) const {
                return (short_options.count(short_opt) > 0);
            }

            const vector<string>& GetAllArguments() const {
                return arguments;
            }
    };

}
