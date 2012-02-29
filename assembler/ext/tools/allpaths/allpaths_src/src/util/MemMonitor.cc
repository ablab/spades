#include "MainTools.h"

#include <fstream>
#include <signal.h>
#include <time.h>

using namespace std;

// convertion factors between jiffies and seconds (check your system!!!)
const double sec_jif = 100.0;
const double jif_sec = 1.0 / sec_jif;


istream &GetLine(istream &is, String &Str, char delimiter = '\n') {
	string str;
	getline(is, str);
	
	Str = str;

	return is;
}

// gets a token from a line of space separated tokens
String LineToken(const String &S, int i) {
  string s = S;
  if (i == 0) return s.substr(0, s.find(" "));
  int j0 = s.find(" (") + 2;
  return s.substr(j0, s.find(")") - j0);
}

// gets a token from a line (from /proc/<pid>/stat) of space separated tokens
String StatLineToken(const String &S, int i) {
  string s = S;
  if (i == 0) return s.substr(0, s.find(" "));
  if (i == 1) {
    int j0 = s.find(" (") + 2;
    return s.substr(j0, s.find(")") - j0);
  }
  i -= 2;
  int j0 = s.find(") ") + 2;
  while (i--) { j0 = s.find(" ", j0) + 1; }
  int j1 = s.find(" ", j0);
  return s.substr(j0, j1 - j0);
}






class Process {
public:
  Process(int pid = 0) : _pid(pid), _cmdline(""), 
                         _etime(0), _etime_prev(0), _detime(0),
                         _utime(0), _utime_prev(0), _dutime(0),
                         _vmsize(0), _vmlck(0), _vmrss(0), 
                         _vmdata(0), _vmstk(0), _vmexe(0), _vmlib(0) {}

  bool IsRunning() { return IsDirectory("/proc/" + ToString(_pid)); }


  void Load() {
    String line;

    ifstream statusstream(("/proc/" + ToString(_pid) + "/status").c_str());
    while (GetLine(statusstream, line)) {
      if (line.Contains("Vm")) {
        String field = line.Before(":");
        int value = WhiteSpaceFree(line.After(":").Before("kB")).Int();

        if      (field == "VmSize") { _vmsize = value; }
        else if (field == "VmLck" ) { _vmlck  = value; }
        else if (field == "VmRSS" ) { _vmrss  = value; }
        else if (field == "VmData") { _vmdata = value; }
        else if (field == "VmStk" ) { _vmstk  = value; }
        else if (field == "VmExe" ) { _vmexe  = value; }
        else if (field == "VmLib" ) { _vmlib  = value; }
      }
    }
    statusstream.close();

    

    ifstream statstream(("/proc/" + ToString(_pid) + "/stat").c_str());
    GetLine(statstream, line);
    statstream.close();

    // utime is the user space time in seconds (i.e. cpu time)

    _utime_prev = _utime;
    _utime = StatLineToken(line, 13).Double() * jif_sec;   // the 13th entry on /proc/<pid>/stat
    _dutime = _utime - _utime_prev;
    

    // etime is the elapsed in seconds

    double t_start_s = StatLineToken(line, 21).Double() * jif_sec;  // the 21st entry on /proc/<pid>/stat

    // get the time since boot
    ifstream uptimestream("/proc/uptime");
    GetLine(uptimestream, line);
    uptimestream.close();
    double t_now_s = LineToken(line, 0).Double();  // the 1st entry on /proc/uptime
    
    _etime_prev = _etime;
    _etime = t_now_s - t_start_s;
    
    // Bug fix: sometimes t_start_s is not properly loaded.
    // In this case, just keep _etime constant, so _detime = 0.
    if ( t_start_s == 0 ) _etime = _etime_prev;
    
    _detime = _etime - _etime_prev;

  }

  String GetTraceback() {
    String traceback;
    String tracefile = "/tmp/traceback_from_process_" + ToString(_pid);
    Remove(tracefile);

    int kill_status = kill(_pid, SIGUSR1);
    while(!IsRegularFile(tracefile)) { usleep(10000); }

    String line;
    ifstream tracestream(tracefile.c_str());
    while(GetLine(tracestream, line)) { traceback += "\n#\t" + line;	}
    tracestream.close();

    Remove(tracefile);
    return traceback;
  }
  
  void Abort() { kill(_pid, SIGABRT); }
 
  String DataLineHeader() {
    return "time(sec), etime, detime, utime, dutime, vmsize(KB), vmrss, vmdata, vmstk, vmexe, vmlib";
  }

  String DataLine() {
    return ToString(time(NULL)) + 
      " " + ToString(_etime) + 
      " " + ToString(_detime) + 
      " " + ToString(_utime) + 
      " " + ToString(_dutime) + 
      " " + ToString(_vmsize) + 
      " " + ToString(_vmrss) + 
      " " + ToString(_vmdata) + 
      " " + ToString(_vmstk) + 
      " " + ToString(_vmexe) + 
      " " + ToString(_vmlib);
  }


  String CmdLine() {
    if (_cmdline == "") {
      ifstream cmdlinestream(("/proc/" + ToString(_pid) + "/cmdline").c_str());

      char cmdchar;
      while (cmdlinestream.get(cmdchar)) {
        _cmdline += (cmdchar == 0) ? ' ' : cmdchar;
      }
      cout << endl;
      
      cmdlinestream.close();
    }
    
    return _cmdline;
  }

  
  double &ETime()  { return _etime;  }
  double &DETime()  { return _detime;  }

  double &UTime()  { return _utime;  }
  double &DUTime()  { return _dutime;  }

  int &VmSize()  { return _vmsize; }
  int &VmLck()   { return _vmlck;  }
  int &VmRss()   { return _vmrss;  }
  int &VmData()  { return _vmdata; }
  int &VmStk()   { return _vmstk;  }
  int &VmExe()   { return _vmexe;  }
  int &VmLib()   { return _vmlib;  }

private:
  int _pid;
  String _cmdline;

  double _etime;           // elapsed time since process started
  double _etime_prev;
  double _detime;

  double _utime;           // time in user mode (calculated from /proc/<pid>/stat[13])
  double _utime_prev;
  double _dutime;

  int _vmsize;
  int _vmlck;
  int _vmrss;
  int _vmdata;
  int _vmstk;
  int _vmexe;
  int _vmlib;


};


int main(int argc, char **argv) {
  RunTimeNoTraceback();

  BeginCommandArguments;
  CommandArgument_Int_Doc(PID,                                  "Process ID of the program to monitor.");
  CommandArgument_Int_OrDefault_Doc(INTERVAL,               20, "The time (in seconds) between output of stats to file.");
  CommandArgument_Int_OrDefault_Doc(MEM_LIMIT,               0, "Memory limit (in kB) of process to watch.  If the process exceeds this limit, it will be killed.");
  CommandArgument_Bool_OrDefault_Doc(WITH_TRACEBACK,         0, "Should a traceback be printed everytime the program's memory usage is polled?");
  CommandArgument_Bool_OrDefault_Doc(TRACEBACK_ON_OVERFLOW,  1, "Should a traceback be printed just prior to killing a program that has exceeded the memory limit?");
  CommandArgument_Bool_OrDefault_Doc(SUMMARY,                0, "Only print the maximum memory usage once the process completes.");
  CommandArgument_String_OrDefault_Doc(OUT,      "/dev/stderr", "Place to write output.  If this is a directory, write the output to 'OUT/PID.mm.log'.");
  EndCommandArguments;

  Process maxps;
  Process ps(PID);

  
  if (ps.IsRunning()) {

    ofstream fout( IsDirectory(OUT) ? (OUT + "/" + ToString(PID) + ".mm.log").c_str() : OUT.c_str() );

    fout << "# (MM," << PID << ") cmdline: '" << ps.CmdLine() << "'" << endl;
    fout << "# (MM," << PID << ") fields: " << ps.DataLineHeader() << endl;

    int iter = 0;
    do {
      ps.Load();
      
      if (iter++ % INTERVAL == 0) {
        fout << "(MM," << PID << ") "  << ps.DataLine() << endl << flush;
        
        if (WITH_TRACEBACK) {
          fout << "# (MM) process traceback: " << ps.GetTraceback() << endl;
        }
      }

      if (MEM_LIMIT > 0 && ps.VmSize() > MEM_LIMIT) {
        fout << "# (MM) notice: monitored process exceeded memory threshold of " << MEM_LIMIT << " kB; killing " << PID << endl;
        
        if (TRACEBACK_ON_OVERFLOW && !WITH_TRACEBACK) {
          fout << "# (MM) process final traceback: " << ps.GetTraceback() << endl;
        }
        
        ps.Abort();
      }
      
      maxps.UTime() = ps.UTime();
      maxps.ETime() = ps.ETime();
      if (ps.DUTime() > maxps.DUTime()) { maxps.DUTime() = ps.DUTime(); }
      if (ps.DETime() > maxps.DETime()) { maxps.DETime() = ps.DETime(); }
      if (ps.VmSize() > maxps.VmSize()) { maxps.VmSize() = ps.VmSize(); }
      if (ps.VmRss()  > maxps.VmRss())  { maxps.VmRss()  = ps.VmRss();  }
      if (ps.VmData() > maxps.VmData()) { maxps.VmData() = ps.VmData(); }
      if (ps.VmStk()  > maxps.VmStk())  { maxps.VmStk()  = ps.VmStk();  }
      if (ps.VmExe()  > maxps.VmExe())  { maxps.VmExe()  = ps.VmExe();  }
      if (ps.VmLib()  > maxps.VmLib())  { maxps.VmLib()  = ps.VmLib();  }
        
      sleep(1);
    } while (ps.IsRunning());

    fout << "# (MM," << PID << ") fields: " << maxps.DataLineHeader() << endl;
    fout << "# (MM," << PID << ") Summary: " << maxps.DataLine() << endl;

    fout.close();
  }

  return 0;
}
