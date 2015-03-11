#include <chrono>
#include <unordered_map>
#include <string>
#include <iostream>
#include <vector>
#include <algorithm>

struct Timers
{

  typedef std::chrono::high_resolution_clock clock;
  typedef clock::duration duration;
  typedef clock::time_point time_point;

  struct Timer
  {

    bool _running;
    duration _duration;
    time_point _start_time;
    int _nCalls;

    Timer()
      : _running(false)
      , _duration(duration::zero())
      , _nCalls(0)
    {}

    void start()
    {
      if (!_running)
        _start_time = clock::now();
      _running = true;
      ++_nCalls;
    }

    bool running() const
    {
      return _running;
    }

    int nCalls() const
    {
      return _nCalls;
    }

    void stop()
    {
      if (_running)
        {
          time_point end_time = clock::now();
          _duration += (end_time - _start_time);
          _running = false;
        }
    }

    duration elapsed() const
    {
      return _duration;
    }

    void reset(bool stop = true)
    {
      if (stop)
        {
          _duration = duration::zero();
          _running = false;
        }
      else
        {
          _duration = duration::zero();
          _running = true;
          _start_time = clock::now();
        }
    }

  };

  void start(std::string name)
  {
    _timers[name].start();
  }

  void stop(std::string name)
  {
    _timers[name].stop();
  }

  bool running(std::string name)
  {
    return _timers[name].running();
  }

  int nCalls(std::string name)
  {
    return _timers[name].nCalls();
  }

  duration elapsed(std::string name) const
  {
    return _timers[name].elapsed();
  }

  static double to_seconds(duration d)
  {
    return std::chrono::duration_cast<std::chrono::duration<double> >(d).count();
  }

  double seconds(std::string name) const
  {
    return to_seconds(_timers[name].elapsed());
  }

  void print_timer(std::ostream& os, std::string name, double scaling_factor = 0) const
  {
    double s = seconds(name);
    os << name << ": " << s;
    if (scaling_factor > 0)
      os << " scaled: " << s/scaling_factor;
    os << std::endl;
  }

  void print_timer_per_call(std::ostream& os, std::string name) const
  {
    print_timer(os, name, _timers[name].nCalls());
  }

  void print_timers(std::ostream& os, double scaling_factor = 0)
  {
    std::vector<std::string> names;
    names.reserve(_timers.size());
    for (auto& e : _timers)
      names.push_back(e.first);
    std::sort(names.begin(),names.end());
    for (auto name : names)
      {
        auto s = seconds(name);
        os << name << ": " << s;
        if (scaling_factor > 0)
          os << " scaled: " << s/scaling_factor;
        os << std::endl;
      }
  }

  void print_timers_per_call(std::ostream& os)
  {
    std::vector<std::string> names;
    names.reserve(_timers.size());
    for (auto& e : _timers)
      names.push_back(e.first);
    std::sort(names.begin(),names.end());
    for (auto name : names)
      {
        int scaling_factor = _timers[name].nCalls();
        auto s = seconds(name);
        os << name << ": " << s;
        if (scaling_factor > 0)
          os << " scaled: " << s/scaling_factor;
        os << std::endl;
      }
  }

  void reset(bool stop = true)
  {
    for (auto& e : _timers)
      e.second.reset(stop);
  }

  mutable std::unordered_map<std::string,Timer> _timers;

};

Timers timers;
