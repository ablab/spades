/*
GetOpt_pp:	Yet another C++ version of getopt.
    Copyright (C) 2007, 2008  Daniel Gutson, FuDePAN
    
    This file is part of GetOpt_pp.

    GetOpt_pp is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    board-games is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef GETOPT_PP_H
#define GETOPT_PP_H

#include <string>
#include <vector>
#include <map>
#include <sstream>

/*
	DESIGN GOALS:
		- EASY to use
		- EASY to learn
		- mimc STL's streams
		- EASY to extend
*/

#ifndef GETOPT_INLINE
#	define GETOPT_INLINE
#endif

namespace GetOpt {

typedef std::vector<std::string> OptionArgs;

struct OptionData
{
	enum _Flags
	{
		CmdLine_NotExtracted,
		CmdLine_Extracted,
		Envir
	};
	
	_Flags flags;
	OptionArgs args;
	OptionData() : flags(CmdLine_NotExtracted) {}
	void clear() { flags = CmdLine_NotExtracted; args.clear(); }
};

typedef std::map<std::string, OptionData> LongOptions;
typedef std::map<char, OptionData> ShortOptions;


struct _Option
{
	enum Result
	{
		OK,
		ParsingError,
		OptionNotFound,
		BadType,
		NoArgs,
		TooManyArgs,
		OptionNotFound_NoEx
	};

	virtual Result operator() (ShortOptions& short_ops, LongOptions& long_ops, std::ios::fmtflags flags) const = 0;
	virtual ~_Option(){}
};

template <class T> inline _Option::Result convert(const std::string& s, T& result, std::ios::fmtflags flags)
{
	std::stringstream ss;
	ss.clear();
	ss.flags(flags);
	ss << s;
	ss >> result;
	if (ss.fail() || !ss.eof())
		return _Option::BadType;
	else
		return _Option::OK;
}

template <> inline _Option::Result convert<std::string>(const std::string& s, std::string& result, std::ios::fmtflags flags)
{
	result = s;
	return _Option::OK;
}


template <class T> class _OptionTBase : public _Option
{
	const char short_opt;
	const std::string long_opt;
protected:
	T& target;
	virtual Result _assign(const OptionArgs& args, std::ios::fmtflags flags) const = 0;
	
public:
	_OptionTBase(const _OptionTBase<T>& other)
		: short_opt(other.short_opt), long_opt(other.long_opt), target(other.target)
	{}
	
	_OptionTBase(char short_opt, const std::string& long_opt, T& target)
		: short_opt(short_opt), long_opt(long_opt), target(target)
	{}
	
	virtual Result operator() (ShortOptions& short_ops, LongOptions& long_ops, std::ios::fmtflags flags) const
	{
		Result ret = OptionNotFound;
		ShortOptions::iterator it = short_ops.find(short_opt);
		
		if (it != short_ops.end())
		{
			it->second.flags = OptionData::CmdLine_Extracted;
			ret = _assign(it->second.args, flags);
		}
		else if (!long_opt.empty())
		{
			LongOptions::iterator it = long_ops.find(long_opt);
			if (it != long_ops.end())
			{
				it->second.flags = OptionData::CmdLine_Extracted;
				ret = _assign(it->second.args, flags);
			}
		}

		return ret;
	}
};


template <class T> class _OptionT : public _OptionTBase<T>
{
protected:
	virtual _Option::Result _assign(const OptionArgs& args, std::ios::fmtflags flags) const
	{		
		switch (args.size())
		{
			case 0:
				return _Option::NoArgs;
				
			case 1:
				return convert<T>(args[0], this->target, flags);

			default:
				return _Option::TooManyArgs;
		}
			
	}
public:	
	_OptionT(const _OptionT<T>& other)
		: _OptionTBase<T>(other)
	{}
	
	_OptionT(char short_opt, const std::string& long_opt, T& target)
		: _OptionTBase<T>(short_opt, long_opt, target)
	{}

};

template <class T> class _OptionT<std::vector<T> > : public _OptionTBase<std::vector<T> >
{
protected:
	virtual _Option::Result _assign(const OptionArgs& args, std::ios::fmtflags flags) const
	{
		if (!args.empty())
		{
			_Option::Result result;
			OptionArgs::const_iterator it = args.begin();
			T temp;
			
			do
			{
				result = convert<T>(*it, temp, flags);
				if (result == _Option::OK)
					this->target.push_back(temp);
					
				++it;
			} while(it != args.end() && result == _Option::OK);
			
			return result;
		}
		else
			return _Option::NoArgs;
	}
	
public:	
	_OptionT(const _OptionT<std::vector<T> >& other)
		: _OptionTBase<std::vector<T> >(other)
	{}

	_OptionT(char short_opt, const std::string& long_opt, std::vector<T>& target)
		: _OptionTBase<std::vector<T> >(short_opt, long_opt, target)
	{}
};


template <class T, class BaseOption>
class _DefValOption : public BaseOption
{
	const T default_value;
public:
	
	_DefValOption(const _DefValOption<T, BaseOption>& other)
		: BaseOption(other), default_value(other.default_value)
	{}

	_DefValOption(char short_opt, const std::string& long_opt, T& target, const T& default_value)
		: BaseOption(short_opt, long_opt, target), default_value(default_value)
	{}

	virtual _Option::Result operator() (ShortOptions& short_ops, LongOptions& long_ops, std::ios::fmtflags flags) const
	{
		_Option::Result ret = BaseOption::operator()(short_ops, long_ops, flags);
		
		if (ret == _Option::OptionNotFound)
		{
			this->target = default_value;
			ret = _Option::OK;
		}

		return ret;
	}
};


template <class T>
inline _OptionT<T> Option(char short_opt, const std::string& long_opt, T& target)
{
	return _OptionT<T>(short_opt, long_opt, target);
}

template <class T>
inline _OptionT<T> Option(char short_opt, T& target)
{
	return _OptionT<T>(short_opt, std::string(), target);
}

// Defaulted version
template <class T>
inline _DefValOption<T, _OptionT<T> > 
Option(char short_opt, const std::string& long_opt, T& target, const T& def)
{
	return _DefValOption<T, _OptionT<T> >(short_opt, long_opt, target, def);
}

template <class T>
inline _DefValOption<T, _OptionT<T> > Option(char short_opt, T& target, const T& def)
{
	return _DefValOption<T, _OptionT<T> >(short_opt, std::string(), target, def);
}

// Defaults for strings:
inline _DefValOption<std::string, _OptionT<std::string> > 
Option(char short_opt, const std::string& long_opt, std::string& target, const char* def)
{
	return _DefValOption<std::string, _OptionT<std::string> >(short_opt, long_opt, target, def);
}

inline _OptionT<std::string> Option(char short_opt, std::string& target, const char* def)
{
	return _DefValOption<std::string, _OptionT<std::string> >(short_opt, std::string(), target, def);
}


class OptionPresent : public _Option
{
	const char short_opt;
	const std::string long_opt;
	bool* const present;
public:
	// two combinations: with/without target, and with/without long opt.
	
	// WITH long_opt:
	OptionPresent(char short_opt, const std::string& long_opt, bool& present)
		: short_opt(short_opt), long_opt(long_opt), present(&present)
	{}

	OptionPresent(char short_opt, const std::string& long_opt)
		: short_opt(short_opt), long_opt(long_opt), present(NULL)
	{}
	
	// WITHOUT long_opt:
	OptionPresent(char short_opt, bool& present)
		: short_opt(short_opt), present(&present)
	{}

	OptionPresent(char short_opt)
		: short_opt(short_opt), present(NULL)
	{}
	
protected:
	virtual Result operator() (ShortOptions& short_ops, LongOptions& long_ops, std::ios::fmtflags flags) const
	{
		bool found;
		ShortOptions::iterator it = short_ops.find(short_opt);
		
		found = (it != short_ops.end());
		if (found)
		{
			it->second.flags = OptionData::CmdLine_Extracted;
		}
		else if (!long_opt.empty())
		{
			LongOptions::iterator it = long_ops.find(long_opt);
			found = (it != long_ops.end());
			if (found)
			{
				it->second.flags = OptionData::CmdLine_Extracted;
			}
		}
		
		if (present != NULL)
		{
			*present = found;
			return OK;
		}
		else
		{
			return found ? OK : OptionNotFound_NoEx;
		}
	}
};

class GetOptEx : public std::exception {};
struct ParsingErrorEx : GetOptEx{};
struct InvalidFormatEx : GetOptEx{};
struct ArgumentNotFoundEx : GetOptEx{};
struct TooManyArgumentsEx : GetOptEx{};
struct OptionNotFoundEx : GetOptEx{};
struct TooManyOptionsEx : GetOptEx{};

enum _EnvTag
{
	Include_Environment
};

class GetOpt_pp
{
	ShortOptions _shortOps;
	LongOptions _longOps;
	std::ios_base::iostate _exc;
	_Option::Result _last;
	std::ios::fmtflags _flags;
	std::string _app_name;

	GETOPT_INLINE void _init_flags();
	GETOPT_INLINE void _parse(int argc, char* argv[]);
	GETOPT_INLINE void _parse_env();
public:
	static const char EMPTY_OPTION;
	
	GETOPT_INLINE GetOpt_pp(int argc, char* argv[]);
	GETOPT_INLINE GetOpt_pp(int argc, char* argv[], _EnvTag);
	
	std::ios_base::iostate exceptions ( ) const			{ return _exc; }
	void exceptions ( std::ios_base::iostate except )	{ _exc = except; }
	void exceptions_all()                               { _exc = std::ios_base::failbit | std::ios_base::eofbit; }
	
	operator bool() const								{ return _last == _Option::OK;	}

	GETOPT_INLINE bool options_remain() const;
	
	void end_of_options() const throw(GetOptEx)
	{
    	if (options_remain() && (_exc & std::ios_base::eofbit))
    	    throw TooManyOptionsEx();
    }
	
	std::ios::fmtflags flags() const					{ return _flags; }
	void flags(std::ios::fmtflags flags)				{ _flags = flags; }
	
	const std::string& app_name() const					{ return _app_name; }
	
	GETOPT_INLINE GetOpt_pp& operator >> (const _Option& opt) throw(GetOptEx);
	
	GETOPT_INLINE GetOpt_pp& operator >> (std::ios_base& (*iomanip)(std::ios_base&));

	// Alternative to manipulators, for those who don't like them: the 'getopt' method :)	
	// Combination 1: with long option:
	template <class T> inline T getopt(char short_opt, const std::string& long_opt) throw(GetOptEx)
	{
		T result;
		operator >> (Option(short_opt, long_opt, result));
		return result;
	}

	template <class T> inline T getopt(char short_opt, const std::string& long_opt, const T& def_value)
	{
		T result;
		operator >> (Option(short_opt, long_opt, result, def_value));
		return result;
	}

	// Combination 2: without long option:
	template <class T> inline T getopt(char short_opt) throw(GetOptEx)
	{
		T result;
		operator >> (Option(short_opt, result));
		return result;
	}

	template <class T> inline T getopt(char short_opt, const T& def_value)
	{
		T result;
		operator >> (Option(short_opt, result, def_value));
		return result;
	}

	typedef std::pair<ShortOptions::const_iterator, LongOptions::const_iterator> ItPair;
	
	template <class Container, class Adapter, class OptionType>
	class _iterator
	{
		typename Container::const_iterator _it;
	public:
		_iterator(ItPair p)
		{
			_it = Adapter::adapt(p);
		}
		
		_iterator(){}
		
		_iterator<Container, Adapter, OptionType>& operator = (const _iterator<Container, Adapter, OptionType>& other)
		{
			_it = other._it;
			return *this;
		}
		
		bool operator != (const _iterator<Container, Adapter, OptionType>& other) const	{	return _it != other._it;	}
		
		OptionType option() const			{	return _it->first;				}
		const OptionArgs&  args() const		{	return _it->second.args;		}
		_iterator<Container, Adapter, OptionType>& operator ++()	{	++_it; return *this;	}
	};
	
	ItPair begin() const
	{
		return ItPair(_shortOps.begin(), _longOps.begin() );
	}
	
	ItPair end() const
	{
		return ItPair(_shortOps.end(), _longOps.end());
	}
	
	struct ShortAdapter
	{
		static ShortOptions::const_iterator adapt(ItPair p)
		{
			return p.first;
		}
	};
	
	struct LongAdapter
	{
		static LongOptions::const_iterator adapt(ItPair p)
		{
			return p.second;
		}
	};
	
	typedef _iterator<ShortOptions, ShortAdapter, char> short_iterator;
	typedef _iterator<LongOptions, LongAdapter, const std::string&> long_iterator;
};

class Environment
{
	// Coming soon!
};

}

#endif
